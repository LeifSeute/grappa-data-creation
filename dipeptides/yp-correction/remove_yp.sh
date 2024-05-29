n which the script is located
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Define the source and destination directories
SOURCE_DIR="$SCRIPT_DIR/../data"
DEST_DIR="$SCRIPT_DIR/../data_removed"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find directories with 'Y' or 'P' in their names and move them to data_removed
find "$SOURCE_DIR" -type d \( -name '*Y*' -o -name '*P*' \) -exec mv {} "$DEST_DIR" \;

