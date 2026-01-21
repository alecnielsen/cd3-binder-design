#!/usr/bin/env bash
#
# Ralph Wiggum Loop: Scientific Code Review with Fresh Context
#
# This script implements iterative code review where each iteration starts
# with fresh context (no memory of previous iterations). This helps catch
# issues that might be overlooked when context accumulates.
#
# Usage:
#   ./ralph_review.sh <module>          # Review a single module
#   ./ralph_review.sh all               # Review all modules (sequential)
#   ./ralph_review.sh all --parallel    # Review all modules in parallel
#   ./ralph_review.sh --status          # Show current tracking status
#
# Environment variables:
#   MAX_ITERATIONS  - Maximum iterations per module (default: 5)
#   CLAUDE_MODEL    - Model to use (default: claude-sonnet-4-20250514)
#   VERBOSE         - Set to 1 for verbose output

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
TRACKING_FILE="$SCRIPT_DIR/tracking.yaml"
PROMPTS_DIR="$SCRIPT_DIR/prompts"
SRC_DIR="$REPO_ROOT/src"

MAX_ITERATIONS="${MAX_ITERATIONS:-5}"
CLAUDE_MODEL="${CLAUDE_MODEL:-claude-sonnet-4-20250514}"
VERBOSE="${VERBOSE:-0}"

# Available modules
# - src/ modules: analysis, formatting, design, structure, pipeline
# - root modules: scripts, modal
MODULES=(analysis formatting design structure pipeline scripts modal)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

    # Check for PyYAML
    if ! python3 -c "import yaml" 2>/dev/null; then
        log_error "PyYAML not installed. Install with: pip install pyyaml"
        exit 1
    fi
}

# Get list of available modules
get_modules() {
    printf '%s\n' "${MODULES[@]}"
}

# Check if module is valid
is_valid_module() {
    local module="$1"
    for m in "${MODULES[@]}"; do
        if [[ "$m" == "$module" ]]; then
            return 0
        fi
    done
    return 1
}

# Get module README path
get_module_readme() {
    local module="$1"
    # scripts and modal are root-level directories
    if [[ "$module" == "scripts" || "$module" == "modal" ]]; then
        echo "$REPO_ROOT/$module/README.md"
    else
        echo "$SRC_DIR/$module/README.md"
    fi
}

# Read tracking YAML value
read_tracking() {
    local key="$1"
    python3 -c "
import yaml
with open('$TRACKING_FILE') as f:
    data = yaml.safe_load(f)
keys = '$key'.split('.')
val = data
for k in keys:
    val = val.get(k, '')
print(val if val else '')
"
}

# Update tracking YAML (with file locking for parallel safety)
update_tracking() {
    local module="$1"
    local field="$2"
    local value="$3"
    local lock_dir="$TRACKING_FILE.lock"

    # Use mkdir for atomic locking (portable, works on macOS)
    local max_attempts=50
    local attempt=0
    while ! mkdir "$lock_dir" 2>/dev/null; do
        ((attempt++))
        if [[ "$attempt" -ge "$max_attempts" ]]; then
            log_warning "Could not acquire lock for tracking.yaml after $max_attempts attempts"
            break
        fi
        sleep 0.1
    done

    # Ensure lock is released on exit
    trap "rmdir '$lock_dir' 2>/dev/null" EXIT

    python3 << EOF
import yaml
from datetime import datetime

with open('$TRACKING_FILE', 'r') as f:
    data = yaml.safe_load(f)

# Update module field
if '$field' == 'iterations':
    data['modules']['$module']['$field'] = int('$value')
elif '$field' == 'issues_found':
    # Append to list
    if not data['modules']['$module'].get('issues_found'):
        data['modules']['$module']['issues_found'] = []
    data['modules']['$module']['issues_found'].append('$value')
else:
    data['modules']['$module']['$field'] = '$value'

# Update timestamp
if '$field' in ['status', 'last_result']:
    data['modules']['$module']['last_reviewed'] = datetime.now().isoformat()

# Update summary
total = 0
clean = 0
issues = 0
for m in data['modules'].values():
    total += m.get('iterations', 0)
    if m.get('status') == 'clean':
        clean += 1
    elif m.get('status') == 'issues':
        issues += 1

data['summary']['total_reviews'] = total
data['summary']['modules_clean'] = clean
data['summary']['modules_with_issues'] = issues
data['summary']['last_run'] = datetime.now().isoformat()

with open('$TRACKING_FILE', 'w') as f:
    yaml.dump(data, f, default_flow_style=False, sort_keys=False)
EOF

    # Release lock
    rmdir "$lock_dir" 2>/dev/null || true
}

# Get source directory for a module
get_module_source_dir() {
    local module="$1"
    # scripts and modal are root-level directories
    if [[ "$module" == "scripts" || "$module" == "modal" ]]; then
        echo "$REPO_ROOT/$module"
    else
        echo "$SRC_DIR/$module"
    fi
}

# Collect source code for a module
collect_source_code() {
    local module="$1"
    local output=""
    local module_dir
    module_dir=$(get_module_source_dir "$module")

    if [[ -d "$module_dir" ]]; then
        while IFS= read -r -d '' pyfile; do
            local rel_path="${pyfile#$REPO_ROOT/}"
            output+="
=== FILE: $rel_path ===
$(cat "$pyfile")

"
        done < <(find "$module_dir" -name "*.py" -type f -print0 | sort -z)
    fi

    # For structure module, also include modal/ files (for context)
    if [[ "$module" == "structure" ]]; then
        local modal_dir="$REPO_ROOT/modal"
        if [[ -d "$modal_dir" ]]; then
            while IFS= read -r -d '' pyfile; do
                local rel_path="${pyfile#$REPO_ROOT/}"
                output+="
=== FILE: $rel_path ===
$(cat "$pyfile")

"
            done < <(find "$modal_dir" -name "*.py" -type f -print0 | sort -z)
        fi
    fi

    echo "$output"
}

# Run a single review iteration
run_review_iteration() {
    local module="$1"
    local iteration="$2"

    log_info "Running review iteration $iteration for module: $module"

    # Get the prompt
    local prompt_file="$PROMPTS_DIR/${module}.md"
    if [[ ! -f "$prompt_file" ]]; then
        log_error "Prompt file not found: $prompt_file"
        return 1
    fi

    # Get module README context
    local readme_path
    readme_path=$(get_module_readme "$module")
    local readme_context=""
    if [[ -f "$readme_path" ]]; then
        readme_context=$(cat "$readme_path")
    else
        log_warning "Module README not found: $readme_path"
    fi

    # Collect source code
    local source_code
    source_code=$(collect_source_code "$module")

    # Read the prompt
    local prompt
    prompt=$(cat "$prompt_file")

    # Check for iteration history log
    local history_file="$SCRIPT_DIR/logs/${module}_history.md"
    local history_context=""
    if [[ -f "$history_file" ]]; then
        history_context=$(cat "$history_file")
    fi

    # Build history section if we have previous iterations
    local history_section=""
    if [[ -n "$history_context" ]]; then
        history_section="
---

# PREVIOUS ITERATION HISTORY

The following shows what previous review iterations found and fixed.
Use this to understand what has already been addressed and avoid repeating the same findings.

$history_context

---
"
    fi

    # Construct the full review request
    local full_request="
$prompt
$history_section
---

# MODULE README CONTEXT

The following is the scientific context from the module's README.md file.
If you find this context is wrong or incomplete, you should update it directly in the README.

$readme_context

---

# SOURCE CODE TO REVIEW

$source_code

---

Please review the source code above against the README context and the scientific criteria in the prompt.
Remember: If no code issues are found, respond with exactly NO_ISSUES (and nothing else).
Note: Updating the module README.md to fix documentation issues does NOT count as a code issue.
"

    # Create a temporary file for the request
    local temp_file
    temp_file=$(mktemp)
    echo "$full_request" > "$temp_file"

    # Run claude with --print for non-interactive mode
    local result
    if [[ "$VERBOSE" == "1" ]]; then
        log_info "Sending request to Claude..."
    fi

    # Use claude --print with full permissions for automated review loop
    result=$(claude --print --dangerously-skip-permissions < "$temp_file" 2>&1) || {
        log_error "Claude review failed"
        rm -f "$temp_file"
        return 1
    }

    rm -f "$temp_file"

    # Log this iteration's result to history file
    local timestamp
    timestamp=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    {
        echo "## Iteration $iteration ($timestamp)"
        echo ""
        echo "$result"
        echo ""
        echo "---"
        echo ""
    } >> "$history_file"

    # Check for NO_ISSUES - robust parsing that handles whitespace and extra output
    # Look for NO_ISSUES on its own line (with possible whitespace)
    if echo "$result" | grep -qE '^\s*NO_ISSUES\s*$'; then
        log_success "Module $module: NO_ISSUES found"
        update_tracking "$module" "status" "clean"
        update_tracking "$module" "last_result" "NO_ISSUES"
        echo ""  # Return empty hash to signal success
        return 0
    else
        log_warning "Module $module: Issues found"
        update_tracking "$module" "status" "issues"
        update_tracking "$module" "last_result" "issues_found"

        # Print the issues to stdout
        echo ""
        echo "=== Issues found in $module (iteration $iteration) ==="
        echo "$result"
        echo "========================================================"
        echo ""

        # Return hash of issues for stuck detection
        echo "$result" | shasum -a 256 | cut -d' ' -f1
        return 1
    fi
}

# Review a single module (loops until clean or max iterations)
review_module() {
    local module="$1"

    log_info "Starting review loop for module: $module"

    # Check if module exists
    if ! is_valid_module "$module"; then
        log_error "Unknown module: $module"
        log_info "Available modules: $(get_modules | tr '\n' ' ')"
        return 1
    fi

    # Get current iteration count
    local current_iterations
    current_iterations=$(read_tracking "modules.${module}.iterations")
    current_iterations=${current_iterations:-0}

    # Check if already at max
    if [[ "$current_iterations" -ge "$MAX_ITERATIONS" ]]; then
        log_warning "Module $module has already reached max iterations ($MAX_ITERATIONS)"
        return 0
    fi

    # Update status to in_progress
    update_tracking "$module" "status" "in_progress"

    # Stuck detection - track previous issue hash
    local prev_issue_hash=""
    local stuck_count=0
    local max_stuck=2  # Exit if same issues seen this many times

    # Loop until clean or max iterations
    while [[ "$current_iterations" -lt "$MAX_ITERATIONS" ]]; do
        local iteration=$((current_iterations + 1))
        update_tracking "$module" "iterations" "$iteration"

        log_info "Module $module: iteration $iteration/$MAX_ITERATIONS"

        # Run review and capture output (including issue hash on failure)
        local output
        local exit_code
        output=$(run_review_iteration "$module" "$iteration")
        exit_code=$?

        if [[ "$exit_code" -eq 0 ]]; then
            # Clean - we're done
            log_success "Module $module passed review (iteration $iteration)"
            return 0
        fi

        # Issues found - check for stuck condition
        local issue_hash
        issue_hash=$(echo "$output" | tail -1)  # Hash is on last line

        if [[ "$issue_hash" == "$prev_issue_hash" ]]; then
            ((stuck_count++))
            log_warning "Module $module: same issues detected ($stuck_count/$max_stuck)"

            if [[ "$stuck_count" -ge "$max_stuck" ]]; then
                log_error "Module $module appears stuck - same issues $max_stuck times in a row"
                log_error "Human intervention may be needed"
                update_tracking "$module" "status" "stuck"
                return 1
            fi
        else
            stuck_count=0
        fi

        prev_issue_hash="$issue_hash"
        current_iterations=$iteration

        # Brief pause between iterations to avoid hammering API
        sleep 2
    done

    log_warning "Module $module reached max iterations ($MAX_ITERATIONS) with issues remaining"
    return 1
}

# Review all modules sequentially
review_all_sequential() {
    log_info "Reviewing all modules sequentially..."

    local failed=0
    for module in "${MODULES[@]}"; do
        if ! review_module "$module"; then
            ((failed++)) || true
        fi
    done

    echo ""
    if [[ "$failed" -eq 0 ]]; then
        log_success "All ${#MODULES[@]} modules passed review!"
    else
        log_warning "Review complete. Modules with issues: $failed"
    fi

    # Return 0 on success (all clean), 1 if any failed
    [[ "$failed" -eq 0 ]]
}

# Review modules in batches (respects Claude subscription concurrency limits)
# Usage: review_all_parallel [batch_size]
# Default batch_size is 2 (safe for most subscriptions)
review_all_parallel() {
    local batch_size="${BATCH_SIZE:-2}"
    log_info "Reviewing all modules in batches of $batch_size..."

    local modules_list=("${MODULES[@]}")
    local temp_dir
    temp_dir=$(mktemp -d)
    local total_failed=0

    # Process in batches
    local i=0
    while [[ "$i" -lt "${#modules_list[@]}" ]]; do
        local batch=()
        local pids=()

        # Build batch
        for ((j=0; j<batch_size && i+j<${#modules_list[@]}; j++)); do
            batch+=("${modules_list[i+j]}")
        done

        log_info "Starting batch: ${batch[*]}"

        # Launch batch in background
        for module in "${batch[@]}"; do
            (
                review_module "$module" > "$temp_dir/${module}.log" 2>&1
                echo $? > "$temp_dir/${module}.exit"
            ) &
            pids+=($!)
        done

        # Wait for batch to complete
        wait "${pids[@]}"

        # Report batch results
        for module in "${batch[@]}"; do
            local exit_code
            exit_code=$(cat "$temp_dir/${module}.exit" 2>/dev/null || echo "1")

            if [[ "$exit_code" == "0" ]]; then
                log_success "$module: clean"
            else
                log_warning "$module: issues found"
                ((total_failed++)) || true
            fi
        done

        ((i+=batch_size))
    done

    # Show detailed logs
    echo ""
    for module in "${modules_list[@]}"; do
        echo "=== Results for $module ==="
        cat "$temp_dir/${module}.log" 2>/dev/null || echo "(no output)"
        echo "==========================="
        echo ""
    done

    # Cleanup
    rm -rf "$temp_dir"

    if [[ "$total_failed" -eq 0 ]]; then
        log_success "All ${#modules_list[@]} modules passed review!"
    else
        log_warning "Review complete. Modules with issues: $total_failed"
    fi

    [[ "$total_failed" -eq 0 ]]
}

# Review all modules (entry point)
review_all() {
    local parallel="${1:-false}"

    if [[ "$parallel" == "true" ]]; then
        review_all_parallel
    else
        review_all_sequential
    fi
}

# Show current status
show_status() {
    echo ""
    echo "=== Review Tracking Status ==="
    echo ""

    python3 - "$TRACKING_FILE" << 'PYEOF'
import sys
import yaml

tracking_file = sys.argv[1]

with open(tracking_file, 'r') as f:
    data = yaml.safe_load(f)

print("Modules:")
for name, info in data.get('modules', {}).items():
    status = info.get('status', 'unknown')
    iterations = info.get('iterations', 0)
    last_reviewed = info.get('last_reviewed', 'never')
    if last_reviewed and last_reviewed != 'null' and last_reviewed is not None:
        last_reviewed = str(last_reviewed)[:19]  # Truncate ISO timestamp
    else:
        last_reviewed = 'never'

    status_emoji = {
        'pending': '[ ]',
        'in_progress': '[~]',
        'clean': '[+]',
        'issues': '[!]'
    }.get(status, '[?]')

    print(f"  {status_emoji} {name}: {status} (iterations: {iterations}, last: {last_reviewed})")

print("")
summary = data.get('summary', {})
print(f"Summary:")
print(f"  Total reviews: {summary.get('total_reviews', 0)}")
print(f"  Modules clean: {summary.get('modules_clean', 0)}")
print(f"  Modules with issues: {summary.get('modules_with_issues', 0)}")
print(f"  Last run: {summary.get('last_run', 'never')}")
PYEOF
}

# Main entry point
main() {
    check_dependencies

    if [[ $# -eq 0 ]]; then
        echo "Usage: $0 <module|all|--status> [--parallel]"
        echo ""
        echo "Modules: $(get_modules | tr '\n' ' ')"
        echo ""
        echo "Options:"
        echo "  <module>     Review a specific module"
        echo "  all          Review all modules (sequential)"
        echo "  all --parallel  Review all modules in parallel"
        echo "  --status     Show current tracking status"
        echo ""
        echo "Environment variables:"
        echo "  MAX_ITERATIONS  Max iterations per module (default: 5)"
        echo "  VERBOSE         Set to 1 for verbose output"
        exit 1
    fi

    case "$1" in
        --status|-s)
            show_status
            ;;
        all)
            # Check for --parallel flag
            if [[ "${2:-}" == "--parallel" ]]; then
                review_all "true"
            else
                review_all "false"
            fi
            ;;
        *)
            review_module "$1"
            ;;
    esac
}

main "$@"
