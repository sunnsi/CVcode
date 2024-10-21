# Abstract
In the CVcode project, we employ the finite difference method to solve partial differential equations, considering both scenarios with and without the MOS constraint. Our code generates scattered data for peak current density ($I_{pc}$), based on varying diffusion coefficients ($D_s$), electrode radius ($R_0$), and MOS parameters ($C_{max}$). We then fit this data to the analytical formula for $I_{pc}$ derived in our article, using the differential evolution method. Finally, we compare the scattered data with the corresponding formula, both with and without the MOS constraint, to validate the accuracy of our results.
### code file
- **spherical_Cmax.py**: Implements the finite difference method for solving partial differential equations under the MOS constraint.
- **spherical_noCmax.py**: Implements the finite difference method for solving partial differential equations without the MOS constraint for comparison.
- **datastorageforMOS.py**: Handles data storage for scenarios under the MOS constraint by function calls.
- **datastoragefornoMOS.py**: Handles data storage for scenarios without the MOS constraint by function calls.
- **co_fitforI.py**: Utilizes the differential evolution method to fit the coefficients for the derived analytical formula of $I_{pc}$.
- **totalfitforMOS.py**: Matches the scattered data of $I_{pc}$ under the MOS constraint with the derived formula.
- **totalfitfornoMOS.py**: Matches the scattered data of $I_{pc}$ without the MOS constraint with the traditional formula from Bard's textbook for comparison.
### data file
- **dataforMOS.xls**: Contains scattered data for $I_{pc}$ under the MOS constraint.
- **datafornoMOS.xls**: Contains scattered data for $I_{pc}$ without the MOS constraint.
- **datalineforMOS.xls**: Contains data generated by the derived formula for $I_{pc}$ under the MOS constraint.
- **datalinefornoMOS.xls**: Contains data generated by the formula for $I_{pc}$ from Bard's electrochemical textbook.
- **dataR2forI.xls**: Compares data from **dataforMOS.xls** and **datalineforMOS.xls** to calculate the $R^2$ index.
