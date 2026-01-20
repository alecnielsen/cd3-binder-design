#!/usr/bin/env bash
#
# Ralph Wiggum Loop: Scientific Code Review with Fresh Context
#
# This script implements iterative code review where each iteration starts
# with fresh context (no memory of previous iterations). This helps catch
# issues that might be overlooked when context accumulates.
#
# Usage:
#   ./run_review.sh <module>          # Review a single module
#   ./run_review.sh all               # Review all modules
#   ./run_review.sh --status          # Show current tracking status
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
SLICES_FILE="$SCRIPT_DIR/context_slices.yaml"
PROMPTS_DIR="$SCRIPT_DIR/prompts"
README_UPDATES="$SCRIPT_DIR/readme_updates.md"
README_PATH="$REPO_ROOT/README.md"

MAX_ITERATIONS="${MAX_ITERATIONS:-5}"
CLAUDE_MODEL="${CLAUDE_MODEL:-claude-sonnet-4-20250514}"
VERBOSE="${VERBOSE:-0}"

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
    python3 "$SCRIPT_DIR/extract_context.py" --list 2>/dev/null | \
        grep -E "^  [a-z]+:" | \
        sed 's/://g' | \
        awk '{print $1}'
}

# Extract README context for a module
extract_context() {
    local module="$1"
    python3 "$SCRIPT_DIR/extract_context.py" "$module"
}

# Get source paths for a module
get_source_paths() {
    local module="$1"
    python3 "$SCRIPT_DIR/extract_context.py" "$module" --sources-only
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

# Update tracking YAML
update_tracking() {
    local module="$1"
    local field="$2"
    local value="$3"

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
}

# Collect source code for a module
collect_source_code() {
    local module="$1"
    local output=""

    while IFS= read -r src_path; do
        local full_path="$REPO_ROOT/$src_path"
        if [[ -d "$full_path" ]]; then
            # It's a directory, collect all .py files
            while IFS= read -r -d '' pyfile; do
                local rel_path="${pyfile#$REPO_ROOT/}"
                output+="
=== FILE: $rel_path ===
$(cat "$pyfile")

"
            done < <(find "$full_path" -name "*.py" -type f -print0 | sort -z)
        elif [[ -f "$full_path" ]]; then
            # It's a file
            local rel_path="${full_path#$REPO_ROOT/}"
            output+="
=== FILE: $rel_path ===
$(cat "$full_path")

"
        fi
    done < <(get_source_paths "$module")

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

    # Extract README context
    local readme_context
    readme_context=$(extract_context "$module")

    # Collect source code
    local source_code
    source_code=$(collect_source_code "$module")

    # Read the prompt
    local prompt
    prompt=$(cat "$prompt_file")

    # Construct the full review request
    local full_request="
$prompt

---

# README CONTEXT (relevant sections)

$readme_context

---

# SOURCE CODE TO REVIEW

$source_code

---

Please review the source code above against the README context and the scientific criteria in the prompt.
Remember: If no issues are found, respond with exactly NO_ISSUES (and nothing else).
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

    # Use claude --print with the input piped
    result=$(claude --print < "$temp_file" 2>&1) || {
        log_error "Claude review failed"
        rm -f "$temp_file"
        return 1
    }

    rm -f "$temp_file"

    # Check for NO_ISSUES
    if echo "$result" | grep -q "^NO_ISSUES$"; then
        log_success "Module $module: NO_ISSUES found"
        update_tracking "$module" "status" "clean"
        update_tracking "$module" "last_result" "NO_ISSUES"
        return 0
    else
        log_warning "Module $module: Issues found"
        update_tracking "$module" "status" "issues"
        update_tracking "$module" "last_result" "issues_found"

        # Append to readme_updates.md
        {
            echo ""
            echo "## Review: $module (iteration $iteration)"
            echo "**Date**: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
            echo ""
            echo "$result"
            echo ""
            echo "---"
        } >> "$README_UPDATES"

        log_info "Issues logged to $README_UPDATES"
        return 1
    fi
}

# Review a single module
review_module() {
    local module="$1"

    log_info "Starting review for module: $module"

    # Check if module exists
    if ! get_modules | grep -q "^${module}$"; then
        log_error "Unknown module: $module"
        log_info "Available modules: $(get_modules | tr '\n' ' ')"
        return 1
    fi

    # Get current iteration count
    local current_iterations
    current_iterations=$(read_tracking "modules.${module}.iterations")
    current_iterations=${current_iterations:-0}

    # Check if we've hit max iterations
    if [[ "$current_iterations" -ge "$MAX_ITERATIONS" ]]; then
        log_warning "Module $module has reached max iterations ($MAX_ITERATIONS)"
        return 0
    fi

    # Update status to in_progress
    update_tracking "$module" "status" "in_progress"

    # Run the review
    local iteration=$((current_iterations + 1))
    update_tracking "$module" "iterations" "$iteration"

    if run_review_iteration "$module" "$iteration"; then
        # Clean - no need for more iterations
        log_success "Module $module passed review (iteration $iteration)"
        return 0
    else
        # Issues found - could run more iterations if under limit
        if [[ "$iteration" -lt "$MAX_ITERATIONS" ]]; then
            log_info "Module $module has issues. Re-run to continue reviewing (iteration $iteration/$MAX_ITERATIONS)"
        else
            log_warning "Module $module reached max iterations with issues remaining"
        fi
        return 1
    fi
}

# Review all modules
review_all() {
    log_info "Reviewing all modules..."

    local modules
    modules=$(get_modules)

    local failed=0
    for module in $modules; do
        if ! review_module "$module"; then
            ((failed++)) || true
        fi
    done

    echo ""
    log_info "Review complete. Modules with issues: $failed"

    return $failed
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
        echo "Usage: $0 <module|all|--status>"
        echo ""
        echo "Modules: $(get_modules | tr '\n' ' ')"
        echo ""
        echo "Options:"
        echo "  <module>   Review a specific module"
        echo "  all        Review all modules"
        echo "  --status   Show current tracking status"
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
            review_all
            ;;
        *)
            review_module "$1"
            ;;
    esac
}

main "$@"
