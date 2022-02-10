package FEM;

public class ModelCode {

	public static void main(String[] args) {
		Node n1 = new Node(1, 0, 0);
		Node n2 = new Node(0, 1, 0);
		Constraint c = new Constraint(true, false, true);
		Force f1 = new Force(1.2, -4, 0);
		Element e1 = new Element(2.1e8, 0.2, n1, n2);

		// set force and constraint
		n1.setConstraint(c);
		n2.setForce(f1);
		// print

		System.out.println("Constraint (u1, u2, u3)");
		n1.getConstraint().print();

		System.out.println("Force (f1, f2, f3)");
		n2.getForce().print();

		System.out.println("Position of node 1 (n11, n12, n13)");
		n1.print();
		System.out.println("Position of node 2 (n21, n22, n23)");
		n2.print();
		System.out.println("Element (E,A,L)");
		e1.print();
	}

}
