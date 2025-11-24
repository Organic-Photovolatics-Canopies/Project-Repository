- **Design Diagrams:**
  - **Level 0: Context Diagram**
    ![Context Diagram](./D0%20-%20Context%20Diagram.png) _Description:_ Shows
    the overall system context, including the user, fabrication process, and
    design query flow. Conventions: UML actor for user, rounded rectangles for
    processes. Purpose: To illustrate the main external interactions with the
    OPV design system.

  - **Level 1: Container Diagram**
    ![Container Diagram](./D1%20-%20Container%20Diagram.png) _Description:_
    Breaks down the system into main containers: Front End (Web App), OPV Design
    Tool, and their interactions. Conventions: UML actor for user, containers as
    rounded rectangles. Purpose: To show how the user interacts with the web app
    and how prompts and responses flow through the system.

  - **Level 2: Component Diagram**
    ![Component Diagram](./D2%20-%20Component%20Diagram.png) _Description:_
    Details the internal components such as Materials List, Graph Neural
    Network, and their data flows (e.g., absorption curves, generated design).
    Conventions: Cylinders for data stores, rounded rectangles for components.
    Purpose: To describe the internal structure and data flow of the OPV design
    model.
