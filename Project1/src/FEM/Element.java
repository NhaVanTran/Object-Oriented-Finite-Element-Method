package FEM;

import inf.jlinalg.Array2DMatrix;
import inf.jlinalg.BLAM;
import inf.jlinalg.IMatrix;
import inf.jlinalg.Vector3D;

public class Element {
	private double area;
	private double eModulus;
	private int[] dofNumbers = new int[6];
	private Node node1;
	private Node node2;

	/* A truss element has two nodes, area, E, dofNumbers */
	public Element(double e, double a, Node n1, Node n2) {
		this.area = a;
		this.eModulus = e;
		this.node1 = n1;
		this.node2 = n2;
	}

	public IMatrix computeStiffnessMatrix() {
		double E = this.eModulus;
		double A = this.area;
		double L = this.getLength();
		// Transformation Matrix
		// T= | getE1 0 0 0 | getE1 = [X2-X1, Y2-Y1, Z2-Z1]/Length(Node1,Node2)
		// ...| 0 0 0 get1E1| getE1 = [m, n, p]
		double m = this.getE1().getX1();
		double n = this.getE1().getX2();
		double p = this.getE1().getX3();
		double[][] transfer = { { m, n, p, 0, 0, 0 }, { 0, 0, 0, m, n, p } };
		IMatrix T = new Array2DMatrix(transfer);
		// Semi-Local stiffness matrix ^^
		double[][] ke_l = { { 1, -1 }, { -1, 1 } };
		IMatrix ke_local = new Array2DMatrix(ke_l);
		IMatrix tmp = new Array2DMatrix(6, 2);
		IMatrix KE_global = new Array2DMatrix(6, 6);
		BLAM.multiply(E * A / L, BLAM.TRANSPOSE, T, BLAM.NO_TRANSPOSE,
				ke_local, 0.0, tmp); // tmp = (E*A/L)*Transpose(T)*ke_local
		BLAM.multiply(1.0, BLAM.NO_TRANSPOSE, tmp, BLAM.NO_TRANSPOSE, T, 0.0,
				KE_global); // KE_global = tmp*T
		return KE_global;
	}

	public Vector computeElementInternalForce(Vector u) {
		// Referent position
		Vector3D X1 = this.node1.getPosition();
		Vector3D X2 = this.node2.getPosition();
		Vector X = new Vector(X1.getX1(), X1.getX2(), X1.getX3(), X2.getX1(),
				X2.getX2(), X2.getX3());
		// Current position
		Vector x = X.add(1, u);

		double E = this.eModulus;
		double A = this.area;
		double L = this.getLength();
		double l = Math.sqrt(Math.pow(x.get(0) - x.get(3), 2)
				+ Math.pow(x.get(1) - x.get(4), 2)
				+ Math.pow(x.get(2) - x.get(5), 2));
		double a = x.get(0) - x.get(3);
		double b = x.get(1) - x.get(4);
		double c = x.get(2) - x.get(5);
		Vector r_int = new Vector(a, b, c, -a, -b, -c);
		double constant = E * A * (Math.pow(l, 2) - Math.pow(L, 2))
				/ (2 * Math.pow(L, 3));
//		System.out.println("Constant : " + constant);
		Vector out = r_int.multiply(constant);
		return out;
	}

	public IMatrix computeTangentElementStiffnessMatrix(Vector u) {
		/* Checked */
		// Referent position
		Vector3D X1 = this.node1.getPosition();
		Vector3D X2 = this.node2.getPosition();
		Vector X = new Vector(X1.getX1(), X1.getX2(), X1.getX3(), X2.getX1(),
				X2.getX2(), X2.getX3());
		// Current position
		Vector x = X.add(1, u);

		double E = this.eModulus;
		double A = this.area;
		double L = this.getLength();
		double l = Math.sqrt(Math.pow(x.get(0) - x.get(3), 2)
				+ Math.pow(x.get(1) - x.get(4), 2)
				+ Math.pow(x.get(2) - x.get(5), 2));
		double constant = E * A * (Math.pow(l, 2) - Math.pow(L, 2))
				/ (2 * Math.pow(L, 3));
		// Tensor AT=|+I3 -I3|
		// __________|-I3 +I3|
		double[][] AT = { { 1, 0, 0, -1, 0, 0 }, { 0, 1, 0, 0, -1, 0 },
				{ 0, 0, 1, 0, 0, -1 }, { -1, 0, 0, 1, 0, 0 },
				{ 0, -1, 0, 0, 1, 0 }, { 0, 0, -1, 0, 0, 1 } };

		// Geometric Stiffness Matrix KG=A*E*(l^2-L^2)*AT/(2*L^3)
		double[][] KG = new double[6][6];
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				KG[i][j] = AT[i][j] * constant;
			}
		}
		// Material Stiffness Matrix KM=(E*A/L^3)*C*C_Transpose
		double[][] KM = new double[6][6];
		// Where C=AT*x; C=[6x1]
		// Compute C
		double a = x.get(0) - x.get(3);
		double b = x.get(1) - x.get(4);
		double c = x.get(2) - x.get(5);
		double[] C = { a, b, c, -a, -b, -c };
		// Compute KM. Therefore, Tangent stiffness Matrix
		double[][] KT = new double[6][6];
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				KM[i][j] = E * A * C[i] * C[j] / Math.pow(L, 3);
				KT[i][j] = KG[i][j] + KM[i][j];// To save for loop
			}
		}
		IMatrix KT_e = new Array2DMatrix(KT);
		return KT_e;

	}

	/* Enumerate DOFs for element */
	public void enumerateDOFs() {
		for (int i = 0; i < 3; i++) {
			/* DOF of Node1 */
			this.dofNumbers[i] = this.node1.getDOFNumbers()[i];
			/* DOF of Node2 */
			this.dofNumbers[i + 3] = this.node2.getDOFNumbers()[i];
		}
	}

	public int[] getDOFNumbers() {
		return this.dofNumbers;
	}

	/* Compute Force */
	public double computeForce() {
		// N=E*A*eps(1,1); where eps(1,1)=(l-L)/L
		Vector3D pn1 = this.getNode1().getPosition();
		Vector3D pn2 = this.getNode2().getPosition();
		Vector3D u1 = this.getNode1().getDisplacement();
		Vector3D u2 = this.getNode2().getDisplacement();
		// Position of nodes after deformation
		Vector3D p1 = pn1.add(u1);
		Vector3D p2 = pn2.add(u2);
		// Length of element after deformation
		double l = p1.subtract(p2).norm2();
		return this.getEModulus() * this.getArea() * (l - this.getLength())
				/ this.getLength();
	}

	public Vector3D getE1() {
		// Get uniaxial vector of the two nodes of element
		// (direction from node1 to node2)
		// Subtract the two nodes of element ( node2 - node1 )
		Vector3D vector = this.getNode2().getPosition()
				.subtract(this.getNode1().getPosition());
		return vector.normalize();// Normalize then return
	}

	public Node getNode1() {
		return this.node1;
	}

	public Node getNode2() {
		return this.node2;
	}

	public double getLength() {
		double dx = this.node1.getPosition().getX1()
				- this.node2.getPosition().getX1();
		double dy = this.node1.getPosition().getX2()
				- this.node2.getPosition().getX2();
		double dz = this.node1.getPosition().getX3()
				- this.node2.getPosition().getX3();
		double length = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2)
				+ Math.pow(dz, 2));
		return length;
	}

	public double getArea() {
		return this.area;
	}

	public double getEModulus() {
		return this.eModulus;
	}

	public void print() {
		System.out.println(this.eModulus + "   " + this.area + "   "
				+ this.getLength());
	}

}
