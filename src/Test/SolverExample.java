package Test;

import FEM.*;
import Models.*;
import inf.jlinalg.GeneralMatrixLSESolver;
import inf.jlinalg.ILSESolver;
import inf.jlinalg.IMatrix;
import inf.jlinalg.MatrixFormat;
import inf.jlinalg.QuadraticMatrixInfo;
import inf.text.ArrayFormat;

public class SolverExample {
	public static void main(String[] args) {
		// Structure struct = SmallTetraeder.createStructure();
		// Structure struct = MediumTetraeder.createStructure();
		Structure struct = SmallTetraeder.createStructure();

		// Enumerate Structure
		int DOFs = struct.enumerateDofs();
		System.out.println("Number of equations:   " + DOFs);

		// Print DOF number of nodes
		System.out.println("Node degree of free dom");
		for (int i = 0; i < struct.getNumberOfNodes(); i++) {
			int[] dofNumbers = struct.getNode(i).getDOFNumbers();
			System.out.println(ArrayFormat.format(dofNumbers));
		}

		// Print Element DOFs
		System.out.println("Element degree of free dom");
		for (int i = 0; i < struct.getNumberOfElement(); i++) {
			int[] dofNumbers = struct.getElement(i).getDOFNumbers();
			System.out.println(ArrayFormat.format(dofNumbers));
		}

		// Size of our matrix
		int neq = DOFs;
		// Create the solver object
		ILSESolver solver = new GeneralMatrixLSESolver();
		// Info object for coefficient matrix
		QuadraticMatrixInfo aInfo = solver.getAInfo();
		// Coefficient Matrix
		IMatrix a = solver.getA();
		// Right hand side
		double[] b = new double[neq];
		// Initialize solver
		aInfo.setSize(neq);
		solver.initialize();
		// Assemble Load Vector
		struct.assembleLoadVector(b);
//		for (int i = 0; i < neq; i++) {
//			System.out.println("Force Vector " + i + ": " + b[i]);
//		}
		// Assemble Stiffness matrix
		struct.assembleStiffnessMatrix(a);

		System.out.println("Solving A x = b");
		System.out.println("Matrix A");
		System.out.println(MatrixFormat.format(a));
		System.out.println("Vector b");
		System.out.println(ArrayFormat.format(b));

		try {
			solver.solve(b);
		} catch (Exception e) {
			System.out.println("Solver failed: " + e.getMessage());
		}
		// Print result
		System.out.println("Solution x");
		System.out.println(ArrayFormat.format(b));
	}
}
