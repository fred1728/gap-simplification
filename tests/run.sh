#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
case "${1:-core}" in
  core) exec gap -q -b tests/core.g ;;
  integration) exec gap -q -b tests/integration.g ;;
  all)
    exec gap -q -b tests/core.g tests/integration.g
    ;;
  *) echo "Usage: $0 [core|integration|all]" >&2; exit 2 ;;
esac
