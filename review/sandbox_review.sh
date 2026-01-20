#!/bin/bash
# Sandboxed Ralph Wiggum review runner
#
# Provides:
# - Filesystem isolation (no access to ~/.ssh, ~/.aws, ~/.config, etc.)
# - Only the repo directory is mounted
# - Runs as non-root user
#
# Network isolation note:
# Full network isolation (only allowing api.anthropic.com) requires Linux with iptables.
# On macOS, Docker provides filesystem isolation but network filtering is limited.
# For paranoid mode on macOS, use a Linux VM or remote Linux host.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
IMAGE_NAME="cd3-review-sandbox"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[SANDBOX]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[SANDBOX]${NC} $1"; }
log_error() { echo -e "${RED}[SANDBOX]${NC} $1"; }

# Check for API key
if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
    log_error "ANTHROPIC_API_KEY environment variable not set"
    exit 1
fi

# Build the container if needed
build_container() {
    log_info "Building sandbox container..."
    docker build -t "$IMAGE_NAME" "$SCRIPT_DIR"
}

# Run review in sandbox
run_sandboxed() {
    local args=("$@")

    log_info "Starting sandboxed review..."
    log_warn "Filesystem isolated: only /repo is accessible"
    log_warn "No access to: ~/.ssh, ~/.aws, ~/.config, ~/, /etc, etc."

    # Only use -it flags if we have a TTY (interactive terminal)
    local interactive_flags=""
    if [[ -t 0 && -t 1 ]]; then
        interactive_flags="-it"
    fi

    docker run --rm $interactive_flags \
        --name "cd3-review-$$" \
        -v "$REPO_DIR:/repo:rw" \
        -e "ANTHROPIC_API_KEY=$ANTHROPIC_API_KEY" \
        -e "HOME=/home/reviewer" \
        --workdir /repo \
        --user reviewer \
        --read-only \
        --tmpfs /tmp:rw,noexec,nosuid \
        --tmpfs /home/reviewer:rw,noexec,nosuid \
        --security-opt no-new-privileges \
        --cap-drop ALL \
        "$IMAGE_NAME" \
        ./review/ralph_review.sh "${args[@]}"
}

# Paranoid mode: Linux-only with network filtering
run_paranoid() {
    local args=("$@")

    if [[ "$(uname)" != "Linux" ]]; then
        log_error "Paranoid mode (network filtering) only works on Linux"
        log_error "On macOS, use a Linux VM or accept filesystem-only isolation"
        exit 1
    fi

    log_info "Starting PARANOID mode (network restricted to api.anthropic.com)..."

    # Create isolated network if it doesn't exist
    docker network inspect cd3-isolated >/dev/null 2>&1 || \
        docker network create --internal cd3-isolated

    # Only use -it flags if we have a TTY
    local interactive_flags=""
    if [[ -t 0 && -t 1 ]]; then
        interactive_flags="-it"
    fi

    # Run with network isolation
    # Note: This requires additional iptables rules on the host to allow
    # only api.anthropic.com. See README for setup instructions.
    docker run --rm $interactive_flags \
        --name "cd3-review-$$" \
        --network cd3-isolated \
        -v "$REPO_DIR:/repo:rw" \
        -e "ANTHROPIC_API_KEY=$ANTHROPIC_API_KEY" \
        -e "HOME=/home/reviewer" \
        --workdir /repo \
        --user reviewer \
        --read-only \
        --tmpfs /tmp:rw,noexec,nosuid \
        --tmpfs /home/reviewer:rw,noexec,nosuid \
        --security-opt no-new-privileges \
        --cap-drop ALL \
        "$IMAGE_NAME" \
        ./review/ralph_review.sh "${args[@]}"
}

# Main
case "${1:-}" in
    --build)
        build_container
        ;;
    --paranoid)
        shift
        build_container
        run_paranoid "$@"
        ;;
    --help)
        echo "Usage: $0 [--build|--paranoid] [ralph_review.sh args...]"
        echo ""
        echo "Options:"
        echo "  --build     Build the sandbox container only"
        echo "  --paranoid  Linux-only: restrict network to api.anthropic.com"
        echo ""
        echo "Examples:"
        echo "  $0 --build                    # Build container"
        echo "  $0 analysis                   # Review analysis module"
        echo "  $0 all --parallel             # Review all modules in parallel"
        echo "  $0 --paranoid all             # Linux: full network isolation"
        ;;
    *)
        build_container
        run_sandboxed "$@"
        ;;
esac
