#!/bin/sh
#
# Author: Cory Quammen
# Note: This file will be renamed to PSFEstimator in PSFEstimator.app/Contents/MacOS/
#

echo "Running PSFEstimator executable."

VPSFO_BUNDLE="`echo "$0" | sed -e 's/\/Contents\/MacOS\/PSFEstimator\ [0-9]*.[0-9]*.[0-9]*//'`"
VPSFO_EXECUTABLE="$VPSFO_BUNDLE/Contents/Resources/bin/PSFEstimator"

export "DYLD_LIBRARY_PATH=$VPSFO_BUNDLE/Contents/Resources/lib/"

export "DYLD_FRAMEWORK_PATH=$VPSFO_BUNDLE/Contents/Frameworks/"

exec "$VPSFO_EXECUTABLE"
