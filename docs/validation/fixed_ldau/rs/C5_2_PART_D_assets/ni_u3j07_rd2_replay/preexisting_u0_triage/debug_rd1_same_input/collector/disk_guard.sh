#!/bin/bash
set -u
RUN="$1"
SOLVER_PID="$2"
FLOOR=1073741824
available_bytes() { df -Pk "$RUN" | awk 'NR==2 {print $4 * 1024}'; }
while kill -0 "$SOLVER_PID" 2>/dev/null; do
  available=$(available_bytes)
  if [ "$available" -lt "$FLOOR" ]; then
    printf '{"event":"disk_floor","solver_pid":%s,"available_bytes":%s,"floor_bytes":%s,"action":"SIGTERM"}\n' "$SOLVER_PID" "$available" "$FLOOR" >> "$RUN/disk_guard_events.jsonl"
    kill -TERM "$SOLVER_PID" 2>/dev/null || true
    exit 0
  fi
  sleep 2
done
