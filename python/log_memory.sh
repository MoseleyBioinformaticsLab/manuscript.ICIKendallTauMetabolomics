LOG_FILE="python_log.txt"
INTERVAL=5 # Log every 60 seconds

while true; do
    echo "--- $(date) ---" >> "$LOG_FILE"
    free -m >> "$LOG_FILE"
    sleep "$INTERVAL"
done
