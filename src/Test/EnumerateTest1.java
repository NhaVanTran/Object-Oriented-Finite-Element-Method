package Test;

import inf.text.ArrayFormat;
import FEM.*;

public class EnumerateTest1 {
	public static void main(String[] args) {
		int eqe = 0;
		Node n1 = new Node(0, 0, 0);
		Node n2 = new Node(1, 1, 0);
		Element e1 = new Element(1, 1, n1, n2);
		Constraint c1 = new Constraint(false, true, true);
		n1.setConstraint(c1);
		eqe = n1.enumerateDOFs(eqe);
		eqe = n2.enumerateDOFs(eqe);
		e1.enumerateDOFs();
		System.out.println("Number of eqe" + eqe);
		System.out.println("Eqe numbers of nodal DOFs");
		System.out.println(ArrayFormat.format(n1.getDOFNumbers()));
		System.out.println(ArrayFormat.format(n2.getDOFNumbers()));
		System.out.println("eqe numbers of element DOFs");
		System.out.println(ArrayFormat.format(e1.getDOFNumbers()));
	}

}
