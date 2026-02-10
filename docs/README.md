# OMID USA Documentation

Welcome to the **OMID USA** (AI-Powered Design of Organic Photovoltaics Canopies for Agrivoltaics) documentation. This project uses Graph Neural Networks (GNNs) to predict properties of Organic Photovoltaic (OPV) materials for transparent solar panel applications in agriculture and architecture.

## 📖 Documentation Guide

### Getting Started
- **[Installation Guide](installation.md)** - Set up your development environment
- **[User Guide](user-guide.md)** - Quick start and common workflows
- **[FAQ](faq.md)** - Frequently asked questions

### Technical Reference
- **[User Manual](user-manual.md)** - Comprehensive technical documentation
- **[Data Sources](data-sources.md)** - Dataset information and references
- **[Troubleshooting](troubleshooting.md)** - Common issues and solutions

## 🚀 Quick Links

### External Resources
- **Demo**: [PCE Predictor on Kaggle](https://www.kaggle.com/code/omrjad/pce-predictor) - Interactive notebook demo
- **HuggingFace Demo**: *(Coming Soon)*
- **Primary Dataset**: [OPV2D Dataset](https://github.com/sunyrain/OPV2D)
- **Repository**: [Project Repository](https://github.com/)

### Project Files
- [Project Description](../Project-Description.md)
- [User Stories](../User%20Stories.md)
- [Task List](../Task%20List.md)
- [Design Diagrams](../Design%20Diagrams/)

## 🎯 Project Overview

OMID USA is a machine learning research project that leverages advanced Graph Neural Networks to predict critical properties of organic photovoltaic materials. The system is designed to accelerate the discovery and optimization of transparent solar panels for:

- **Agrivoltaics**: Vineyard and greenhouse canopies that generate power while allowing light transmission for crop growth
- **Building Integration**: Windows and architectural features that combine aesthetics with energy generation
- **Lightweight Applications**: Rooftop installations that don't require structural reinforcement

### Key Features

- **Multi-Task Prediction**: Simultaneously predicts three critical OPV properties:
  - **PCE** (Power Conversion Efficiency) - Energy conversion performance
  - **Voc** (Open-Circuit Voltage) - Maximum voltage output
  - **Jsc** (Short-Circuit Current Density) - Current generation capability

- **Advanced Graph Encoding**: Converts molecular structures (SMILES) into graph representations with sophisticated node and edge features

- **State-of-the-Art Architecture**: Implements Graphormer, a transformer-based GNN architecture for molecular property prediction

- **Comprehensive Pipeline**: End-to-end preprocessing workflow from raw quantum chemistry data to model-ready graph datasets

## 👥 Team Members

- [Dhruv Pratap Singh](../Member%20Bios/Dhruv%20Pratap%20Singh.md) - LLM Integration
- [Ido Gal](../Member%20Bios/Ido%20Gal.md) - GNN Development & Data Processing
- [Milo Ginn](../Member%20Bios/Milo%20Ginn.md) - Front-End Development
- [Om Rajesh Jadhav](../Member%20Bios/Om%20Rajesh%20Jadhav.md) - GNN Development
- [Toan Nham](../Member%20Bios/Toan%20Nham.md) - GenAI & Model Training

## 📊 Current Status

### Completed ✅
- Data acquisition and analysis (OPV2D dataset)
- Complete 10-step preprocessing pipeline
- Baseline GCN implementation (R² = 0.4335, RMSE = 0.7259)
- Graphormer encoding implementation
- Literature review on GNN architectures

### In Progress 🔄
- Graphormer model training and optimization
- TD-DFT absorption data generation for transparency optimization
- Web interface development
- LLM integration for design recommendations
- Model deployment to HuggingFace/Kaggle

## 🔍 What's Next?

1. **Try the Demo**: Check out the [PCE Predictor on Kaggle](https://www.kaggle.com/code/omrjad/pce-predictor) - no installation needed!
2. **New Users**: Start with the [Installation Guide](installation.md) to set up your environment
3. **Quick Start**: Follow the [User Guide](user-guide.md) to run your first preprocessing pipeline
4. **Deep Dive**: Explore the [User Manual](user-manual.md) for detailed technical information
5. **Issues**: Check the [FAQ](faq.md) and [Troubleshooting](troubleshooting.md) guides

> **📝 Note**: Preprocessing notebooks for OPV2D dataset and documentation images will be added soon. 
> See [TODO.md](TODO.md) for planned updates.

## 📝 License & Citation

This project uses data from the OPV2D dataset. If you use this project in your research, please cite appropriately. See [Data Sources](data-sources.md) for citation information.

## 🤝 Contributing

For questions, issues, or contributions, please refer to the main repository or contact the team members listed above.

---

**Last Updated**: February 2026
