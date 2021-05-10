# laminate_analysis
First and Last Ply Failure Load of Laminate using Partial and Full Degradation

## Lamina Class:
    Stores geometric, elastic, strength, hygrothermal properties of the lamina.
    Have functionality to return transformation, compliance, stiffness and reduced stiffness matrix of lamina.

## Laminate Class:
    Stores array of Lamina class object, position of these lamina from mid surface and their failure status.
    Have functionality to return ABD Matrix of the laminate and failure status of laminate.

## Output:
    First and Last Ply failure loads.
    local stresses and local strains in each lamina.
    global stresses and global strains in each lamina.
    Mid surface strains and curvature.