#!/bin/sh
#
# Author: Cory Quammen
# Note: This file will be renamed to VisualPSFOptimizer in VisualPSFOptimizer.app/Contents/MacOS/
#

echo "Running VisualPSFOptimizer executable."

VPSFO_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/VisualPSFOptimizer\ [0-9]*.[0-9]*.[0-9]*//'`"
VPSFO_EXECUTABLE="$QVIGA_BUNDLE/Contents/Resources/bin/VisualPSFOptimizer"

export "DYLD_LIBRARY_PATH=$VPSFO_BUNDLE/Contents/Resources/lib/"

export "DYLD_FRAMEWORK_PATH=$VPSFO_BUNDLE/Contents/Frameworks/"

exec "$VPSFO_EXECUTABLE"
