#!/usr/bin/env bash
#
# Ralph Wiggum Loop: Holistic Scientific Code Review
#
# This script implements iterative holistic code review where each iteration:
# - Starts with fresh context (no conversation memory)
# - Sees the entire codebase structure
# - Uses structured history to track what's been reviewed/fixed
# - Can explore any file using agent-style Read tools
#
# Usage:
#   ./ralph_review.sh              # Run holistic review
#   ./ralph_review.sh --status     # Show current tracking status
#   ./ralph_review.sh --reset      # Reset history and start fresh
#
# Environment variables:
#   MAX_ITERATIONS  - Maximum iterations (default: 10)
#   CLAUDE_MODEL    - Model to use (default: claude-sonnet-4-20250514)
#   CONTEXT_BUDGET  - Max context as fraction (default: 0.5 = 50%)
#   VERBOSE         - Set to 1 for verbose output

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
TRACKING_FILE="$SCRIPT_DIR/tracking.yaml"
PROMPT_FILE="$SCRIPT_DIR/prompts/holistic.md"
HISTORY_FILE="$SCRIPT_DIR/logs/holistic_history.md"
LOGS_DIR="$SCRIPT_DIR/logs"
ITERATION_STATS_FILE="$LOGS_DIR/iteration_stats.tsv"

MAX_ITERATIONS="${MAX_ITERATIONS:-10}"
CLAUDE_MODEL="${CLAUDE_MODEL:-claude-sonnet-4-20250514}"
CONTEXT_BUDGET="${CONTEXT_BUDGET:-0.5}"
VERBOSE="${VERBOSE:-0}"

# Approximate token limits (conservative estimates)
# Claude's context is ~200k tokens; we target CONTEXT_BUDGET of that
MAX_CONTEXT_TOKENS=200000
TARGET_TOKENS=$(python3 -c "print(int($MAX_CONTEXT_TOKENS * $CONTEXT_BUDGET))")

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

log_debug() {
    if [[ "$VERBOSE" == "1" ]]; then
        echo -e "${CYAN}[DEBUG]${NC} $1"
    fi
}

# Check dependencies
check_dependencies() {
    if ! command -v claude &> /dev/null; then
        log_error "claude CLI not found. Install with: npm install -g @anthropic-ai/claude-code"
        exit 1
    fi

    if ! command -v python3 &> /dev/null; then
        log_error "python3 not found"
        exit 1
    fi

    if ! python3 -c "import yaml" 2>/dev/null; then
        log_error "PyYAML not installed. Install with: pip install pyyaml"
        exit 1
    fi
}

# Ensure logs directory exists
ensure_logs_dir() {
    mkdir -p "$LOGS_DIR"
}

# Initialize iteration stats file
init_iteration_stats() {
    echo -e "iteration\tinput_tokens\toutput_tokens\tcontext_pct\tcost_usd\tsummary" > "$ITERATION_STATS_FILE"
}

# Record stats for an iteration
record_iteration_stats() {
    local iteration="$1"
    local input_tokens="$2"
    local output_tokens="$3"
    local cost_usd="$4"
    local summary="$5"

    local context_pct=$(python3 -c "print(f'{100 * $input_tokens / $MAX_CONTEXT_TOKENS:.1f}')")

    echo -e "${iteration}\t${input_tokens}\t${output_tokens}\t${context_pct}%\t\$${cost_usd}\t${summary}" >> "$ITERATION_STATS_FILE"
}

# Extract brief summary from Claude's response
extract_summary() {
    local result="$1"

    # Check for NO_ISSUES
    if echo "$result" | grep -qE '^\s*NO_ISSUES\s*$'; then
        echo "Clean - no issues found"
        return
    fi

    # Count issues fixed
    local issue_count
    issue_count=$(echo "$result" | grep -c "^### Issue" || echo "0")

    if [[ "$issue_count" -gt 0 ]]; then
        # Extract first issue title as example
        local first_issue
        first_issue=$(echo "$result" | grep "^### Issue" | head -1 | sed 's/^### Issue [0-9]*: //')
        echo "${issue_count} issue(s): ${first_issue:0:40}..."
    else
        echo "Review in progress"
    fi
}

# Print summary table at end of loop
print_summary_table() {
    echo ""
    echo "╔════════════════════════════════════════════════════════════════════════════════════╗"
    echo "║                              RALPH LOOP SUMMARY                                    ║"
    echo "╠════════════════════════════════════════════════════════════════════════════════════╣"

    if [[ ! -f "$ITERATION_STATS_FILE" ]]; then
        echo "║  No iteration data available                                                      ║"
        echo "╚════════════════════════════════════════════════════════════════════════════════════╝"
        return
    fi

    printf "║ %-4s │ %-10s │ %-10s │ %-8s │ %-8s │ %-30s ║\n" "Iter" "Input" "Output" "Context" "Cost" "Summary"
    echo "╟──────┼────────────┼────────────┼──────────┼──────────┼────────────────────────────────╢"

    # Skip header line, read data
    tail -n +2 "$ITERATION_STATS_FILE" | while IFS=$'\t' read -r iter input output pct cost summary; do
        # Truncate summary to fit
        summary="${summary:0:30}"
        printf "║ %-4s │ %10s │ %10s │ %8s │ %8s │ %-30s ║\n" "$iter" "$input" "$output" "$pct" "$cost" "$summary"
    done

    echo "╚════════════════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

# Estimate tokens from text (rough: 1 token ≈ 4 characters)
estimate_tokens() {
    local text="$1"
    local chars
    chars=$(echo -n "$text" | wc -c)
    echo $((chars / 4))
}

# Collect ALL source code for review
collect_all_source_code() {
    log_debug "Collecting all source code..."

    local output=""
    local file_count=0
    local total_lines=0

    # Collect all Python files from src/, scripts/, modal/
    for dir in "$REPO_ROOT/src" "$REPO_ROOT/scripts" "$REPO_ROOT/modal"; do
        if [[ -d "$dir" ]]; then
            while IFS= read -r -d '' pyfile; do
                local rel_path="${pyfile#$REPO_ROOT/}"
                local lines
                lines=$(wc -l < "$pyfile" | tr -d ' ')
                ((total_lines += lines))
                ((file_count++))

                output+="
================================================================================
FILE: $rel_path ($lines lines)
================================================================================
$(cat "$pyfile")

"
            done < <(find "$dir" -name "*.py" -type f -print0 2>/dev/null | sort -z)
        fi
    done

    # Also include config.yaml
    if [[ -f "$REPO_ROOT/config.yaml" ]]; then
        local lines
        lines=$(wc -l < "$REPO_ROOT/config.yaml" | tr -d ' ')
        ((total_lines += lines))
        ((file_count++))
        output+="
================================================================================
FILE: config.yaml ($lines lines)
================================================================================
$(cat "$REPO_ROOT/config.yaml")

"
    fi

    log_info "Collected $file_count files, $total_lines total lines"
    echo "$output"
}

# Read tracking value
read_tracking() {
    local key="$1"
    python3 -c "
import yaml
try:
    with open('$TRACKING_FILE') as f:
        data = yaml.safe_load(f) or {}
    keys = '$key'.split('.')
    val = data
    for k in keys:
        val = val.get(k, '')
    print(val if val else '')
except:
    print('')
"
}

# Update tracking file
update_tracking() {
    local field="$1"
    local value="$2"

    python3 << EOF
import yaml
from datetime import datetime

try:
    with open('$TRACKING_FILE', 'r') as f:
        data = yaml.safe_load(f) or {}
except FileNotFoundError:
    data = {}

# Ensure structure exists
if 'holistic' not in data:
    data['holistic'] = {
        'status': 'pending',
        'iterations': 0,
        'last_reviewed': None,
        'files_examined': [],
        'issues_fixed': 0
    }

# Update field
if '$field' == 'iterations':
    data['holistic']['$field'] = int('$value')
elif '$field' == 'issues_fixed':
    data['holistic']['$field'] = int('$value')
else:
    data['holistic']['$field'] = '$value'

# Update timestamp
data['holistic']['last_reviewed'] = datetime.now().isoformat()

with open('$TRACKING_FILE', 'w') as f:
    yaml.dump(data, f, default_flow_style=False, sort_keys=False)
EOF
}

# Initialize or reset tracking
init_tracking() {
    python3 << EOF
import yaml
from datetime import datetime

data = {
    'holistic': {
        'status': 'pending',
        'iterations': 0,
        'last_reviewed': None,
        'files_examined': [],
        'issues_fixed': 0
    },
    'summary': {
        'total_iterations': 0,
        'last_run': datetime.now().isoformat()
    }
}

with open('$TRACKING_FILE', 'w') as f:
    yaml.dump(data, f, default_flow_style=False, sort_keys=False)
EOF
    log_info "Tracking file initialized"
}

# Parse structured history from Claude's output
append_to_history() {
    local iteration="$1"
    local result="$2"
    local timestamp
    timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    {
        echo "## Iteration $iteration ($timestamp)"
        echo ""
        echo "$result"
        echo ""
        echo "---"
        echo ""
    } >> "$HISTORY_FILE"
}

# Get summary of previous iterations for context
get_history_summary() {
    if [[ ! -f "$HISTORY_FILE" ]]; then
        echo "No previous iterations."
        return
    fi

    # Return full history if small enough, otherwise summarize
    local history_size
    history_size=$(wc -c < "$HISTORY_FILE" | tr -d ' ')
    local max_history_chars=$((TARGET_TOKENS * 2))  # ~50% of budget for history

    if [[ "$history_size" -lt "$max_history_chars" ]]; then
        cat "$HISTORY_FILE"
    else
        log_warning "History file large ($history_size bytes), truncating to recent iterations"
        # Get last 3 iterations
        grep -n "^## Iteration" "$HISTORY_FILE" | tail -3 | head -1 | cut -d: -f1 | xargs -I{} tail -n +{} "$HISTORY_FILE"
    fi
}

# Run a single review iteration
run_review_iteration() {
    local iteration="$1"

    log_info "Running holistic review iteration $iteration"

    # Check prompt file exists
    if [[ ! -f "$PROMPT_FILE" ]]; then
        log_error "Prompt file not found: $PROMPT_FILE"
        return 1
    fi

    # Build the review request
    local prompt
    prompt=$(cat "$PROMPT_FILE")

    local source_code
    source_code=$(collect_all_source_code)

    local history
    history=$(get_history_summary)

    # Construct full request
    local full_request="
$prompt

---

# COMPLETE SOURCE CODE

ALL source code is included below. You MUST review every file.

$source_code

---

# PREVIOUS ITERATION HISTORY

Review what has been examined and fixed in previous iterations.
Focus on areas NOT YET REVIEWED or verify that previous fixes are correct.

$history

---

# ITERATION $iteration INSTRUCTIONS

1. Read the history above to understand what's been done
2. Identify areas not yet thoroughly reviewed
3. Use Read tool to examine files - do NOT ask for files to be provided
4. Fix any issues you find by editing files directly
5. Report your findings in the structured format specified in the prompt

Remember: Output NO_ISSUES only if you made zero code edits AND verified previous fixes.
"

    # Estimate context usage
    local prompt_tokens
    prompt_tokens=$(estimate_tokens "$full_request")
    log_debug "Estimated prompt tokens: $prompt_tokens (budget: $TARGET_TOKENS)"

    if [[ "$prompt_tokens" -gt "$TARGET_TOKENS" ]]; then
        log_warning "Prompt exceeds context budget ($prompt_tokens > $TARGET_TOKENS tokens)"
    fi

    # Create temp file for request
    local temp_file
    temp_file=$(mktemp)
    echo "$full_request" > "$temp_file"

    # Run claude with JSON output to get actual token usage
    log_info "Sending request to Claude..."

    local json_output
    if ! json_output=$(claude --print --dangerously-skip-permissions --output-format json < "$temp_file" 2>&1); then
        log_error "Claude review failed"
        rm -f "$temp_file"
        return 1
    fi

    rm -f "$temp_file"

    # Extract result text and usage stats from JSON using Python
    local result
    result=$(echo "$json_output" | python3 -c "import sys, json; print(json.load(sys.stdin).get('result', ''))")

    local input_tokens output_tokens cost_usd
    input_tokens=$(echo "$json_output" | python3 -c "import sys, json; u=json.load(sys.stdin).get('usage',{}); print(u.get('input_tokens',0) + u.get('cache_read_input_tokens',0))")
    output_tokens=$(echo "$json_output" | python3 -c "import sys, json; print(json.load(sys.stdin).get('usage',{}).get('output_tokens',0))")
    cost_usd=$(echo "$json_output" | python3 -c "import sys, json; print(json.load(sys.stdin).get('total_cost_usd',0))")

    local summary
    summary=$(extract_summary "$result")
    record_iteration_stats "$iteration" "$input_tokens" "$output_tokens" "$cost_usd" "$summary"

    # Append to history
    append_to_history "$iteration" "$result"

    # Check for NO_ISSUES
    if echo "$result" | grep -qE '^\s*NO_ISSUES\s*$'; then
        # Verify no uncommitted changes in CODE (exclude review tracking files)
        local code_changes
        code_changes=$(git -C "$REPO_ROOT" status --porcelain -- src/ scripts/ modal/ config.yaml CLAUDE.md README.md)
        if [[ -n "$code_changes" ]]; then
            log_warning "Claude said NO_ISSUES but has uncommitted code changes"
            log_warning "Claude should commit fixes before reporting NO_ISSUES"
            echo "$code_changes"
            echo "$result" | shasum -a 256 | cut -d' ' -f1
            return 1
        fi
        log_success "Holistic review: NO_ISSUES found (verified: code unchanged)"
        update_tracking "status" "clean"
        return 0
    else
        log_warning "Holistic review: Issues found or work done"
        update_tracking "status" "in_progress"

        # Print results
        echo ""
        echo "=== Iteration $iteration Results ==="
        echo "$result"
        echo "===================================="
        echo ""

        # Return hash for stuck detection
        echo "$result" | shasum -a 256 | cut -d' ' -f1
        return 1
    fi
}

# Main review loop
run_holistic_review() {
    log_info "Starting holistic review loop (max $MAX_ITERATIONS iterations)"
    log_info "Context budget: ${CONTEXT_BUDGET} ($TARGET_TOKENS tokens)"

    ensure_logs_dir
    init_iteration_stats

    # Get current iteration count
    local current_iterations
    current_iterations=$(read_tracking "holistic.iterations")
    current_iterations=${current_iterations:-0}

    if [[ "$current_iterations" -ge "$MAX_ITERATIONS" ]]; then
        log_warning "Already at max iterations ($MAX_ITERATIONS). Use --reset to start fresh."
        return 0
    fi

    # Stuck detection
    local prev_issue_hash=""
    local stuck_count=0
    local max_stuck=2

    # Review loop
    while [[ "$current_iterations" -lt "$MAX_ITERATIONS" ]]; do
        local iteration=$((current_iterations + 1))
        update_tracking "iterations" "$iteration"

        log_info "Iteration $iteration/$MAX_ITERATIONS"

        # Run review
        local output
        local exit_code
        output=$(run_review_iteration "$iteration")
        exit_code=$?

        if [[ "$exit_code" -eq 0 ]]; then
            log_success "Review complete - no issues found (iteration $iteration)"
            update_tracking "status" "clean"
            print_summary_table
            return 0
        fi

        # Stuck detection
        local issue_hash
        issue_hash=$(echo "$output" | tail -1)

        if [[ "$issue_hash" == "$prev_issue_hash" ]]; then
            ((stuck_count++))
            log_warning "Same findings detected ($stuck_count/$max_stuck)"

            if [[ "$stuck_count" -ge "$max_stuck" ]]; then
                log_error "Review appears stuck - same findings $max_stuck times"
                log_error "Human intervention may be needed"
                update_tracking "status" "stuck"
                print_summary_table
                return 1
            fi
        else
            stuck_count=0
        fi

        prev_issue_hash="$issue_hash"
        current_iterations=$iteration

        # Pause between iterations
        sleep 2
    done

    log_warning "Reached max iterations ($MAX_ITERATIONS) - review may be incomplete"
    update_tracking "status" "max_iterations"
    print_summary_table
    return 1
}

# Show status
show_status() {
    echo ""
    echo "=== Holistic Review Status ==="
    echo ""

    if [[ ! -f "$TRACKING_FILE" ]]; then
        echo "No tracking file found. Run a review first."
        return
    fi

    python3 - "$TRACKING_FILE" << 'PYEOF'
import sys
import yaml

with open(sys.argv[1], 'r') as f:
    data = yaml.safe_load(f) or {}

holistic = data.get('holistic', {})
status = holistic.get('status', 'unknown')
iterations = holistic.get('iterations', 0)
last_reviewed = holistic.get('last_reviewed', 'never')

status_emoji = {
    'pending': '[ ]',
    'in_progress': '[~]',
    'clean': '[✓]',
    'stuck': '[!]',
    'max_iterations': '[M]'
}.get(status, '[?]')

print(f"Status: {status_emoji} {status}")
print(f"Iterations: {iterations}")
print(f"Last reviewed: {last_reviewed}")

if 'summary' in data:
    print(f"\nTotal iterations (all time): {data['summary'].get('total_iterations', 0)}")
    print(f"Last run: {data['summary'].get('last_run', 'never')}")
PYEOF

    # Show history file size
    if [[ -f "$HISTORY_FILE" ]]; then
        local history_lines
        history_lines=$(wc -l < "$HISTORY_FILE" | tr -d ' ')
        echo ""
        echo "History: $HISTORY_FILE ($history_lines lines)"
    fi
}

# Reset for fresh start
reset_review() {
    log_warning "Resetting holistic review state..."

    # Archive old history if exists
    if [[ -f "$HISTORY_FILE" ]]; then
        local archive_name="$LOGS_DIR/holistic_history_$(date +%Y%m%d_%H%M%S).md"
        mv "$HISTORY_FILE" "$archive_name"
        log_info "Archived old history to: $archive_name"
    fi

    # Reset tracking
    init_tracking

    log_success "Review state reset. Ready for fresh holistic review."
}

# Main entry point
main() {
    check_dependencies

    case "${1:-}" in
        --status|-s)
            show_status
            ;;
        --reset|-r)
            reset_review
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Run holistic scientific code review using Ralph Wiggum loop."
            echo ""
            echo "Options:"
            echo "  (no args)    Run holistic review"
            echo "  --status     Show current review status"
            echo "  --reset      Reset history and start fresh"
            echo "  --help       Show this help"
            echo ""
            echo "Environment variables:"
            echo "  MAX_ITERATIONS   Max iterations (default: 10)"
            echo "  CONTEXT_BUDGET   Context fraction (default: 0.5)"
            echo "  VERBOSE          Set to 1 for debug output"
            ;;
        "")
            run_holistic_review
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
}

main "$@"
