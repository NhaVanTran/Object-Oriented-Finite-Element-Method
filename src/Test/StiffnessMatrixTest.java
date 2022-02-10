package Test;

import inf.jlinalg.IMatrix;
import inf.jlinalg.MatrixFormat;
import FEM.*;

public class StiffnessMatrixTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Node n1 = new Node(0, 0, 0);
		Node n2 = new Node(1, 1, 1);
		Element e = new Element(3, Math.sqrt(3), n1, n2);
		IMatrix ke = e.computeStiffnessMatrix();
		System.out.println("Element Stiffness Matrix");
		System.out.println(MatrixFormat.format(ke));
		Vector u = new Vector(1, 2, 0, 3.22, 2.21, 1.19);
		IMatrix kt = e.computeTangentElementStiffnessMatrix(u);
		System.out.println("Tangention Element Stiffness Matrix");
		System.out.println(MatrixFormat.format(kt));
		Vector r_int = e.computeElementInternalForce(u);
		r_int.print("Internal force = ");

	}

}
