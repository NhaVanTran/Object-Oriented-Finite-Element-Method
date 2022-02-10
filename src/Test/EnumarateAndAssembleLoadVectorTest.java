package Test;

import inf.text.ArrayFormat;
import FEM.*;
import Models.*;

public class EnumarateAndAssembleLoadVectorTest {

	public static void main(String[] args) {
		Structure struct = SmallTetraeder.createStructure();
		// Number of equation
		int DOFs = struct.enumerateDofs();
		System.out.println("Number of equations:   " + DOFs);
		// Print DOF number of nodes
		System.out.println("Node degree of free dom");
		for (int i = 0; i < struct.getNumberOfNodes(); i++) {
			int[] dofNumbers = struct.getNode(i).getDOFNumbers();
			System.out.println(ArrayFormat.format(dofNumbers));
		}
		// Print DOF number of elements
		System.out.println("Element degree of free dom");
		for (int i = 0; i < struct.getNumberOfElement(); i++) {
			int[] dofNumbers = struct.getElement(i).getDOFNumbers();
			System.out.println(ArrayFormat.format(dofNumbers));
		}
		// Test Load vector
		double[] rGlobal = new double[DOFs];
		struct.assembleLoadVector(rGlobal);
		for (int i = 0; i < rGlobal.length; i++) {
			System.out.println("Force Dofs " + i + ": " + rGlobal[i]);

		}
	}

}
