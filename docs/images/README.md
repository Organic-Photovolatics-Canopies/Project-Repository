# Documentation Images

This directory contains images, screenshots, and diagrams for the OMID USA documentation.

## Contents

Place images here that are referenced in the documentation:

- **Screenshots**: Terminal output, GUI interfaces, web demo
- **Diagrams**: Workflow diagrams, architecture visualizations
- **Plots**: Training curves, prediction scatter plots, distribution plots
- **Molecular Visualizations**: Example molecular structures and graphs

## Usage

Reference images in markdown documentation using relative paths:

```markdown
![Pipeline Output](images/pipeline_output.png)
```

## Naming Conventions

Use descriptive, lowercase names with hyphens:

- ✅ `pipeline-output-example.png`
- ✅ `molecular-graph-structure.png`
- ✅ `training-loss-curve.png`
- ❌ `Screenshot 2026-02-10.png`
- ❌ `IMG_1234.png`

## Image Guidelines

- **Format**: PNG for screenshots, SVG for diagrams (when possible)
- **Size**: Optimize images for web (< 500 KB when possible)
- **Resolution**: High enough to read text (at least 1920px wide for screenshots)
- **Alt Text**: Always provide descriptive alt text in markdown

## Examples to Add

Consider adding these images to enhance the documentation:

1. **installation.md**:
   - Successful conda environment creation
   - RDKit version check output
   - Pipeline test run completion

2. **user-guide.md**:
   - Pipeline progress output
   - Output directory structure
   - Molecular graph visualization
   - Training loss curves

3. **user-manual.md**:
   - Architecture diagrams (can reuse from Design Diagrams/)
   - Data flow diagrams
   - Feature engineering visualization

4. **troubleshooting.md**:
   - Common error messages
   - Before/after fixing issues

5. **data-sources.md**:
   - Property distribution plots
   - Dataset statistics visualizations
   - Example molecules from HCEPDB

---

**Note**: Add actual images as the project develops and demo/interface becomes available.
