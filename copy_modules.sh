#!/bin/bash

# Create directories
mkdir -p src/{search,feature_analysis,genetic_diversity,visualization}
mkdir -p tests/{search,feature_analysis,genetic_diversity,visualization,integration}

# Copy source files
cp -v /home/ubuntu/yak-new/src/search/*.py src/search/
cp -v /home/ubuntu/yak-new/src/feature_analysis/*.py src/feature_analysis/
cp -v /home/ubuntu/yak-new/src/genetic_diversity/*.py src/genetic_diversity/
cp -v /home/ubuntu/yak-new/src/visualization/*.py src/visualization/

# Copy test files
cp -v /home/ubuntu/yak-new/tests/search/*.py tests/search/
cp -v /home/ubuntu/yak-new/tests/feature_analysis/*.py tests/feature_analysis/
cp -v /home/ubuntu/yak-new/tests/genetic_diversity/*.py tests/genetic_diversity/
cp -v /home/ubuntu/yak-new/tests/visualization/*.py tests/visualization/
cp -v /home/ubuntu/yak-new/tests/integration/*.py tests/integration/
