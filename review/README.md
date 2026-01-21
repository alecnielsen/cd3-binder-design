# Ralph Wiggum Review System

Automated scientific code review using iterative Claude sessions with fresh context.

## Concept

The "Ralph Wiggum Loop" (named after the Simpsons character) implements iterative AI review where:
- Each iteration starts with **fresh context** (no memory of previous attempts)
- **Files persist** on disk between iterations
- **History accumulates** in a log file, injected into each iteration's prompt
- The iteration that **fixes code cannot exit** - requires verification by fresh session

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     RALPH LOOP                               │
├─────────────────────────────────────────────────────────────┤
│  ┌──────────┐    ┌──────────┐    ┌──────────┐              │
│  │  Spawn   │───▶│  Inject  │───▶│  Attempt │              │
│  │ (fresh)  │    │(prompt+  │    │ (review, │              │
│  └──────────┘    │ code+    │    │  fix,    │              │
│       ▲          │ history) │    │  commit) │              │
│       │          └──────────┘    └────┬─────┘              │
│       │                               │                     │
│  ┌────┴─────┐    ┌──────────┐    ┌────▼─────┐              │
│  │   Loop   │◀───│   Kill   │◀───│ Validate │              │
│  │(if fail) │    │ (wipe    │    │(NO_ISSUES│              │
│  └──────────┘    │ context) │    │+ git ok) │              │
│                  └──────────┘    └──────────┘              │
└─────────────────────────────────────────────────────────────┘
```

## Usage

```bash
# Run the review loop
./ralph_review.sh

# Check current status
./ralph_review.sh --status

# Reset and start fresh (archives old history)
./ralph_review.sh --reset

# Customize
MAX_ITERATIONS=5 ./ralph_review.sh
```

## Files

```
review/
├── ralph_review.sh          # Main loop script
├── prompts/
│   └── holistic.md          # Review criteria and instructions
├── logs/
│   ├── holistic_history.md  # Iteration history (persists across runs)
│   └── iteration_stats.tsv  # Token usage per iteration
└── tracking.yaml            # State (iterations, status)
```

## Exit Conditions

| Condition | Exit Code | Meaning |
|-----------|-----------|---------|
| `NO_ISSUES` + git clean | 0 | Clean verification pass |
| Max iterations reached | 1 | Time limit, may need more work |
| Same findings 2x in a row | 1 | Stuck, needs human intervention |

## Validation Logic

The loop only accepts `NO_ISSUES` if:
1. Claude outputs the literal string `NO_ISSUES`
2. `git status` shows no uncommitted changes in code directories

This ensures: **the iteration that fixes code cannot exit the loop**.
A fresh iteration must verify the fixes before the loop can complete.

## Prompt Structure

Each iteration receives:
1. **Review criteria** (`prompts/holistic.md`) - Scientific requirements
2. **Full source code** - ALL Python files from src/, scripts/, modal/
3. **Iteration history** - What previous iterations found and fixed

## Known Issues

### Context Budget
- Full source code: ~94k tokens
- System prompt: ~23k tokens
- Total input: ~117k tokens (58% of 200k context)
- Limited room for long responses

**Potential fixes:**
- Split repo into smaller review chunks
- Compress code (remove comments, docstrings for review)
- Use agent-style (manifest + Read tool) but enforce coverage

### Token Tracking
The `--output-format json` parsing has edge cases where multiline
results cause parsing failures. Currently catches errors but may
report incorrect token counts.

### Coverage Enforcement
The prompt instructs Claude to review all files and provide a checklist,
but there's no programmatic verification that all files were actually examined.

**Potential fixes:**
- Parse Claude's response for file list, compare against manifest
- Require explicit `[x] filename.py` for each file
- Reject NO_ISSUES if checklist incomplete

## History Format

Each iteration appends to `logs/holistic_history.md`:

```markdown
## Iteration N (2024-01-21T10:00:00Z)

## Issues Found and Fixed
### Issue 1: Description
- **File**: path/to/file.py:LINE
- **Problem**: What was wrong
- **Fix**: What was changed

NO_ISSUES (or issue list)

---
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `MAX_ITERATIONS` | 10 | Maximum loop iterations |
| `CONTEXT_BUDGET` | 0.5 | Target context usage (not enforced) |
| `VERBOSE` | 0 | Set to 1 for debug output |

## References

- [Original Ralph Wiggum technique](https://ghuntley.com/ralph/)
- [Claude Code Ralph plugin](https://github.com/anthropics/claude-code/blob/main/plugins/ralph-wiggum/README.md)
