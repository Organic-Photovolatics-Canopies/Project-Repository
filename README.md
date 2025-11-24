# AI-Powered OPV Design

## Fall Design Report

1. Team Members
   - **Advisor:** Dr. Mohsen Rezayat
   - Milo Ginn
   - Om Jadhav
   - Dhruv Singh
   - Ido Gal
   - Toan Nham

2. Project Description (Assignment #2) See:
   [Project-Description.md](./Project-Description.md) A summary of the project
   scope, goals, and team members.

3. User Stories and Design Diagrams (Assignment #4)
   - **User Stories:** [User_Stories.md](./User_Stories.md) Detailed user
     stories for all stakeholders.

   - **Design Diagrams:**
     - **Level 0: Context Diagram**
       ![Context Diagram](./Design_Diagrams/D0%20-%20Context%20Diagram.png)
       _Description:_ Shows the overall system context, including the user,
       fabrication process, and design query flow. Conventions: UML actor for
       user, rounded rectangles for processes. Purpose: To illustrate the main
       external interactions with the OPV design system.

     - **Level 1: Container Diagram**
       ![Container Diagram](./Design_Diagrams/D1%20-%20Container%20Diagram.png)
       _Description:_ Breaks down the system into main containers: Front End
       (Web App), OPV Design Tool, and their interactions. Conventions: UML
       actor for user, containers as rounded rectangles. Purpose: To show how
       the user interacts with the web app and how prompts and responses flow
       through the system.

     - **Level 2: Component Diagram**
       ![Component Diagram](./Design_Diagrams/D2%20-%20Component%20Diagram.png)
       _Description:_ Details the internal components such as Materials List,
       Graph Neural Network, and their data flows (e.g., absorption curves,
       generated design). Conventions: Cylinders for data stores, rounded
       rectangles for components. Purpose: To describe the internal structure
       and data flow of the OPV design model.
4. Project Tasks and Timeline
5. ABET Concerns
6. [PPT Presentation](./Fall%20Design%20Presentation.pptx)
7. Self-Assessment Essays
   - [Milo Ginn](./Self%20Assessment%20Essays/Milo%20Ginn.md)
8. Professional Biographies
   - [Dhruv Pratap Singh](./Member%20Bios/Dhruv%20Pratap%20Signh.md)
   - [Milo Ginn](./Member%20Bios/Milo%20Ginn.md)
   - [Om Rajesh Jadhav](./Member%20Bios/Om%20Rajesh%20Jadhav.md)
   - [Ido Gal](./Member%20Bios/Ido%20Gal.md)
   - [Toan Nham](./Member%20Bios/Toan%20Nham.md)
9. Appendix

This section provides references, citations, links to code repositories, meeting
notes, and evidence of work for each team member.

**Team Meetings:**

- The CS team met twice weekly for a total of 2.5 hours per week, focusing on
  technical development, data analysis, and code review.
- The full project team also met weekly to coordinate interdisciplinary tasks
  and project milestones.

**Evidence of Effort (45+ hours per member):**

Each member will add their own contributions individually below.

- **Ido Gal:**
  - Work in progress is documented in the `progress/` folder, including:
    - Quantum chemistry dataset analysis
      ([Quantum_Chemistry_Analysis_Report.md](progress/Quantum_Chemistry_Analysis_Report.md))
    - Data cleaning and pipeline documentation
      ([dataset_documentation.md](progress/dataset_documentation.md),
      [pipeline.sql](progress/pipeline.sql))
    - Data files and analysis images (e.g., `data_calcqcset1.csv`,
      `comprehensive_analysis.png`, `database_structure_analysis.png`)
  - Researched the Harvard Clean Energy dataset, performed data cleaning, and
    prepared the data for machine learning model training. This involved:
    - Downloading, exploring, and understanding the dataset structure
    - Cleaning quantum chemistry tables using statistical filtering (see
      `dataset_documentation.md`)
    - Writing SQL and Python scripts to process and validate the data (see
      `pipeline.sql`)
    - Documenting the process and results in progress reports

**References and Citations:**

- Harvard Clean Energy Project Dataset:
  [https://www.census.gov/data/datasets/clean-energy.html](https://www.census.gov/data/datasets/clean-energy.html)
- Quantum Chemistry Data Analysis: See
  `progress/Quantum_Chemistry_Analysis_Report.md`
- Data Cleaning Pipeline: See `progress/dataset_documentation.md` and
  `progress/pipeline.sql`
- Code Repository:
  [Project GitHub](https://github.com/Organic-Photovolatics-Canopies/Project-Repository)
