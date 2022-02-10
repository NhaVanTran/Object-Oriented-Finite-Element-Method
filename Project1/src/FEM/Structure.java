package FEM;

import inf.jlinalg.Array2DMatrix;
import inf.jlinalg.GeneralMatrixLSESolver;
import inf.jlinalg.ILSESolver;
import inf.jlinalg.IMatrix;
import inf.jlinalg.MatrixFormat;
import inf.jlinalg.QuadraticMatrixInfo;
import inf.text.ArrayFormat;

import java.util.ArrayList;

public class Structure {
	private ArrayList<Node> nodes = new ArrayList<Node>();
	private ArrayList<Element> Elements = new ArrayList<Element>();
	private double[] u;// Displacement of System equations

	public Node addNode(double x1, double x2, double x3) {
		Node n = new Node(x1, x2, x3);
		nodes.add(n);
		return n;
	}

	public Element addElement(double e, double a, int n1, int n2) {
		Node nn1 = nodes.get(n1);
		Node nn2 = nodes.get(n2);
		Element ele = new Element(e, a, nn1, nn2);
		Elements.add(ele);
		return ele;
	}

	public int getNumberOfNodes() {
		return nodes.size();
	}

	public Node getNode(int id) {
		return nodes.get(id);
	}

	public int getNumberOfElement() {
		return Elements.size();
	}

	public Element getElement(int id) {
		Element ele = Elements.get(id);
		return ele;
	}

	public int getIndexOfNode(Node n) {
		int ind = 0;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (n == this.getNode(i)) {
				ind = i;
			}
			// break;
		}
		return ind;
	}

	public int enumerateDofs() {
		int nnode = this.nodes.size(); // number of nodes
		// Enumerate for each node
		int[][] DOFs = new int[nnode][3];
		/* DOFs is a matrix [nnode]x[3] */
		int start = 0;
		for (int i = 0; i < nnode; i++) {
			/*
			 * Enumerate DOF number for node then assign that node DOF to total
			 * DOFs
			 */
			start = this.getNode(i).enumerateDOFs(start);
			DOFs[i] = this.getNode(i).getDOFNumbers();
		}
		// Enumerate for each Element
		for (int i = 0; i < this.getNumberOfElement(); i++) {
			this.getElement(i).enumerateDOFs();
		}
		return start;
	}

	// solve
	public void solveLinear() {
		// Enumerate Structure
		int DOFs = this.enumerateDofs();

		// Size of our matrix
		int neq = DOFs;
		// Create the solver object
		ILSESolver solver = new GeneralMatrixLSESolver();
		// Info object for coefficient matrix
		QuadraticMatrixInfo KInfo = solver.getAInfo();
		// Coefficient Matrix
		IMatrix KGlobal = solver.getA();
		// Right hand side
		double[] GGlobal = new double[neq];
		// Initialize solver
		KInfo.setSize(neq);
		solver.initialize();
		// Assemble Load Vector
		this.assembleLoadVector(GGlobal);
		// Assemble Stiffness matrix
		this.assembleStiffnessMatrix(KGlobal);

		// System.out.println("Solving KGlobal u = GGlobal");
		try {
			solver.solve(GGlobal);
		} catch (Exception e) {
			System.out.println("Solver failed: " + e.getMessage());
		}
		// Print result
		this.u = GGlobal;
		System.out.println("Solution of Displacement u");
		System.out.println(ArrayFormat.format(GGlobal));
		this.selectDisplacement(GGlobal);
	}

	public void assembleStiffnessMatrix(IMatrix kGlobal) {
		int nele = this.getNumberOfElement();
		for (int i = 0; i < nele; i++) {
			Element ele = this.getElement(i);
			int[] eleDof = ele.getDOFNumbers();
			IMatrix ke_local = ele.computeStiffnessMatrix();
			for (int j = 0; j < eleDof.length; j++) {
				if (eleDof[j] >= 0) { // Assemble node with Dof >= 0
					for (int j2 = 0; j2 < eleDof.length; j2++) {
						if (eleDof[j2] >= 0) {
							double value = ke_local.get(j, j2)
									+ kGlobal.get(eleDof[j], eleDof[j2]); // Ke[j,j2]+KGlobal[eleDof[j],eleDof[j2]]
							kGlobal.set(eleDof[j], eleDof[j2], value);
						}

					}

				}

			}
		}
	}

	public void assembleLoadVector(double[] rGlobal) {
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			Node node_i = this.getNode(i);
			if (node_i.getForce() != null) { // check if there is force applied
				int[] DOFs = node_i.getDOFNumbers(); // DOF numbers of node_i
				if (node_i.getConstraint() == null) { // No constraint
					for (int j = 0; j < 3; j++) {
						// rGlobal[DOf[j]] = force[i]
						rGlobal[DOFs[j]] = node_i.getForce().getComponent(j);
					}
				} else { // Constrained
					for (int j = 0; j < 3; j++) {
						if (node_i.getConstraint().isFree(j) == true) {
							rGlobal[DOFs[j]] = node_i.getForce()
									.getComponent(j);
						}
					}
				}
			}
		}
	}

	public void selectDisplacement(double[] uGlobal) {
		// Loop of all Nodes
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			Node node = this.getNode(i);
			// DOF of Node. For example [5,6,-1]
			int[] DOFNode = node.getDOFNumbers();
			// Displacement of Node. Example [-0.5,-1,0]
			double[] disNode = new double[3];

			for (int j = 0; j < 3; j++) {
				if (DOFNode[j] != -1) {
					disNode[j] = uGlobal[DOFNode[j]];
				}
			}
			node.setDisplacement(disNode);
		}

	}

	public void printStructure() {
		System.out.println("Number of Elements of The Structure:   "
				+ this.getNumberOfElement());
		System.out.println("Number of Node of The Structure:   "
				+ this.getNumberOfNodes());
	}

	public void printResult() {
		// Print Nodes
		System.out.println("\nListing Structure\n");
		System.out.println("\nNodes");
		System.out.println("idx" + "\t\t\t\t" + "x1" + "\t\t\t\t" + "x2"
				+ "\t\t\t\t" + "x3");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			double[] pnode = this.getNode(i).getPosition().toArray();
			System.out.println(i + "\t\t\t\t" + pnode[0] + "\t\t\t\t"
					+ pnode[1] + "\t\t\t\t" + pnode[2]);
		}

		// Print constraints
		System.out.println("\nConstraints");
		System.out.println("node" + "\t\t\t\t" + "u1" + "\t\t\t\t" + "u2"
				+ "\t\t\t\t" + "u3");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getConstraint() != null) {
				String[] st = new String[3];
				int[] DOF = this.getNode(i).getDOFNumbers();
				for (int j = 0; j < 3; j++) {
					if (DOF[j] != -1) {
						st[j] = "free";
					} else {
						st[j] = "fixed";// DOF[j]==-1 means fixed
					}
				}
				System.out.println(i + "\t\t\t\t" + st[0] + "\t\t\t\t" + st[1]
						+ "\t\t\t\t" + st[2]);
			}

		}

		// Print force
		System.out.println("\nForces");
		System.out.println("node" + "\t\t\t\t" + "r1" + "\t\t\t\t" + "r2"
				+ "\t\t\t\t" + "r3");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getForce() != null) {
				double r1 = this.getNode(i).getForce().getComponent(0);
				double r2 = this.getNode(i).getForce().getComponent(1);
				double r3 = this.getNode(i).getForce().getComponent(2);
				System.out.println(i + "\t\t\t\t" + r1 + "\t\t\t\t" + r2
						+ "\t\t\t\t" + r3);
			}
		}
		// Print Elements
		System.out.println("\nElements");
		System.out.println("idx" + "\t\t\t\t" + "E" + "\t\t\t\t" + "A"
				+ "\t\t\t\t" + "Length");
		for (int i = 0; i < this.getNumberOfElement(); i++) {
			System.out.println(i + "\t\t\t\t"
					+ this.getElement(i).getEModulus() + "\t\t\t\t"
					+ this.getElement(i).getArea() + "\t\t\t\t"
					+ this.getElement(i).getLength());
		}
		// Print analysis result
		System.out.println("\nListing analysis result");
		System.out.println("\nDisplacements");
		System.out.println("node" + "\t\t\t\t" + "u1" + "\t\t\t\t" + "u2"
				+ "\t\t\t\t" + "u3");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			System.out.println(i + "\t\t\t\t"
					+ this.getNode(i).getDisplacement().get(0) + "\t\t\t\t"
					+ this.getNode(i).getDisplacement().get(1) + "\t\t\t\t"
					+ this.getNode(i).getDisplacement().get(2));
		}
		// Print Element Force
		System.out.println("\nElement Force");
		System.out.println("ele" + "\t\t\t\t" + "force");
		for (int i = 0; i < this.getNumberOfElement(); i++) {
			System.out.println(i + "\t\t\t\t"
					+ this.getElement(i).computeForce());
		}
	}

	/***************** Nonlinear Part ********************************/
	public void solveArclengthControl(int step) {
		/*
		 * Using ArcLength Control to solve Non-Linear Problem General sphere
		 * (Crisfield,1981)
		 */
		int DOFs = this.enumerateDofs();
		// Radius of sphere to control
		double s = 1;
		// External force
		double[] r_external = new double[DOFs];
		this.assembleLoadVector(r_external);
		Vector r_ext = new Vector(r_external);
		/* Initial Point */
		Vector u_n = new Vector(DOFs);
		double lamda_n = 0;
		int n = 0;
		while (n < step) {
			/***************** Predictor *****************/
			System.out
					.println("\n\n\n/***************** Predictor *****************/");
			System.out.println("Predictor Step n =    " + n);
			Vector u_k = new Vector(DOFs);
			u_k.assign(u_n);
			double lamda_k = lamda_n;
			// Compute Delta u_Lamda
			Vector Du_lamda = this.computeDuLamda(u_n, r_ext);
			double s0 = Math.sqrt(Du_lamda.scalarProduct(Du_lamda) + 1);
			Vector Du = Du_lamda.multiply(s / s0);
			int sign = 1;
			double checksign = r_ext.scalarProduct(Du) / Du.scalarProduct(Du);
			if (checksign < 0) {
				sign = -1;
			}
			System.out.println("Sign : " + sign);
			Du = Du.multiply(sign);// ///////////
			System.out.println("Norm of Du: " + Du.Norm2());
			double DLamda = sign * s / s0;
			// Update value of displacement and Lamda
			u_k.assign(u_k.add(1, Du));
			lamda_k = lamda_k + DLamda;
			/***************** Corrector iteration (step k) *****************/
			System.out
					.println("/***************** Corrector iteration (step k) *****************/");
			double stopCondition = 1;
			double eps = 1e-5;
			int k = 0;
			while (k < 40 && stopCondition > eps) {
				System.out.println("Step k =    " + k);
				Du_lamda = this.computeDuLamda(u_k, r_ext);
				Vector Dur = this.computeDur(u_k, r_ext, lamda_k);
				DLamda = this.computeDeltaLamda(u_n, u_k, lamda_n, lamda_k, s,
						Du_lamda, Dur);
				Du.assign(Dur.add(DLamda, Dur));
				// Update value of displacement and Lamda
				u_k.assign(u_k.add(1, Du));
				lamda_k = lamda_k + DLamda;
				// Check Convergence
				Vector difference = u_n.add(-1, u_k);
				stopCondition = Du.Norm2() / difference.Norm2();
				k = k + 1;
				System.out.println("Convergence :    " + stopCondition);
				//this.selectDisplacement(u_k.toArray());
			}
			u_n.assign(u_k);
			lamda_n = lamda_k;
			n = n + 1;
			System.out.println("Lamda at step_n: :   " + lamda_n);
			this.u = u_n.toArray();
			//this.selectDisplacement(u_n.toArray());
		}
		this.selectDisplacement(u_n.toArray());
	}

	public double computeDeltaLamda(Vector displacement_n,
			Vector displacement_k, double lamda_n, double lamda_k, double s,
			Vector Du_lamda, Vector Dur) {
		// Copy input to avoid errors related to pointer (^!^)
		Vector u_n = new Vector(displacement_n.getN());
		u_n.assign(displacement_n);
		Vector u_k = new Vector(displacement_k.getN());
		u_k.assign(displacement_k);
		Vector Du_l = new Vector(Du_lamda.getN());
		Du_l.assign(Du_lamda);
		Vector Du_r = new Vector(Dur.getN());
		Du_r.assign(Dur);
		double l_k = lamda_k;
		double l_n = lamda_n;
		// /////////////////////////////////////////////////////////////
		Vector difference = u_k.add(-1, u_n);
		double sk = Math.sqrt(difference.scalarProduct(difference)
				+ Math.pow(l_k - l_n, 2));
		Vector DFu = new Vector(u_n.getN());
		for (int i = 0; i < DFu.getN(); i++) {
			DFu.set(i, (u_k.get(i) - u_n.get(i)) / sk);
		}
		double DFlamda = (l_k - l_n) / sk;
		double Fu = Math.sqrt(difference.scalarProduct(difference)
				+ Math.pow(l_n - l_k, 2))
				- s;
		double numerator = -Fu - DFu.scalarProduct(Du_r);
		double dinomirator = DFu.scalarProduct(Du_l) + DFlamda;

		return numerator / dinomirator;
	}

	public Vector computeDur(Vector displacement, Vector r_external,
			double lamda_k) {
		// Copy input to avoid errors related to pointer (^!^)
		Vector u = new Vector(displacement.getN());
		Vector r_ext = new Vector(displacement.getN());
		u.assign(displacement);
		r_ext.assign(r_external);
		/***********/
		int DOFs = displacement.getN();
		Vector r_int = new Vector(DOFs);
		this.assembleInternalForce(u, r_int);
		IMatrix KT_global = new Array2DMatrix(DOFs, DOFs);
		this.assembleTangentStiffnessMatrix(KT_global, u);
		Vector G = r_int.add(-1 * lamda_k, r_ext);
		Vector rhs = G.multiply(-1);
		return this.solveLinearSystemEquation(KT_global, rhs);
	}

	public void assembleInternalForce(Vector displacement, Vector r_internal) {
		/* This method return the internal force of our structure */
		// Copy input to avoid errors related to pointer (^!^)
		Vector u = new Vector(displacement.getN());
		Vector r_int = new Vector(displacement.getN());
		u.assign(displacement);
		// Copy displacement
		for (int i = 0; i < this.getNumberOfElement(); i++) {
			Element ele = this.getElement(i);
			int[] eleDof = ele.getDOFNumbers();
			// Detect Displacement of this element
			Vector u_ele = new Vector(6);
			for (int j = 0; j < 6; j++) {
				if (eleDof[j] != -1) {
					// Assign displacement values of DOF != -1
					u_ele.set(j, u.get(eleDof[j]));
				}
			}
			// Compute Element vector of internal force
			Vector re = ele.computeElementInternalForce(u_ele);
			// Assembling
			for (int j = 0; j < re.getN(); j++) {
				if (eleDof[j] != -1) {
					// rGlobal[eleDof[j]] = rGlobal[eleDof[j]] + re[j];
					double value = r_int.get(eleDof[j]) + re.get(j);
					r_int.set(eleDof[j], value);
				}
			}

		}

	}

	public Vector computeDuLamda(Vector displacement, Vector r_external) {
		// Copy input to avoid errors related to pointer (^!^)
		Vector u = new Vector(displacement.getN());
		Vector r_ext = new Vector(displacement.getN());
		u.assign(displacement);
		r_ext.assign(r_external);
		/***********/
		int DOFs = displacement.getN();
		IMatrix KT_global = new Array2DMatrix(DOFs, DOFs);
		this.assembleTangentStiffnessMatrix(KT_global, u);
		return this.solveLinearSystemEquation(KT_global, r_ext);
	}

	public Vector solveLinearSystemEquation(IMatrix Kt, Vector r) {
		// Finding the root of system equation Kt * x = r
		// size of our matrix
		int neq = r.getN();
		// Copy right hand side
		Vector x = new Vector(neq);
		x.assign(r); // x = r ;
		double[] root = x.toArray();
		// Create the solver object
		ILSESolver solver = new GeneralMatrixLSESolver();
		// Info object for coefficient matrix
		QuadraticMatrixInfo AInfo = solver.getAInfo();
		// Coefficient Matrix
		IMatrix A = solver.getA();
		// initialize solver
		AInfo.setSize(neq);
		solver.initialize();
		// Set matrix A
		for (int i = 0; i < neq; i++) {
			for (int j = 0; j < neq; j++) {
				A.set(i, j, Kt.get(i, j));
			}
		}
		try {
			solver.solve(root);
		} catch (Exception e) {
			System.out.println("Solver failed: " + e.getMessage());
		}
		return new Vector(root);
	}

	public void assembleTangentStiffnessMatrix(IMatrix KT_global,
			Vector displacement) {
		// Copy displacement
		Vector u = new Vector(displacement.getN());
		u.assign(displacement);
		for (int i = 0; i < this.getNumberOfElement(); i++) {
			Element ele = this.getElement(i);
			int[] eleDof = ele.getDOFNumbers();
			// Detect Displacement of this element
			Vector u_ele = new Vector(6);
			for (int j = 0; j < 6; j++) {
				if (eleDof[j] != -1) {
					// Assign displacement values of DOF != -1
					u_ele.set(j, u.get(eleDof[j]));
				}
			}
			// Compute Tangent Element Stiffness Matrix of element
			IMatrix Kt_local = ele.computeTangentElementStiffnessMatrix(u_ele);
			// Assemble local tangent stiffness matrix to global stiffness
			// matrix
			for (int j = 0; j < eleDof.length; j++) {
				if (eleDof[j] >= 0) { // Assemble node with Dof >= 0
					for (int j2 = 0; j2 < eleDof.length; j2++) {
						if (eleDof[j2] >= 0) {
							double value = Kt_local.get(j, j2)
									+ KT_global.get(eleDof[j], eleDof[j2]); // Kt[j,j2]+KT[eleDof[j],eleDof[j2]]
							KT_global.set(eleDof[j], eleDof[j2], value);
						}

					}

				}

			}
		}
	}

	public void solveNewtonRaphson() {
		int DOFs = this.enumerateDofs();
		double[] rGlobal = new double[DOFs];
		double convergence = 1;
		double eps = 10E-5;
		this.assembleLoadVector(rGlobal);
		Vector r_ext = new Vector(rGlobal);
		double[] u_n = new double[DOFs];
		Vector u_n_vector = new Vector(u_n);
		while (convergence > eps) {
			double[] r_internal = new double[DOFs];
			Vector r_int = new Vector(r_internal);

			u_n = u_n_vector.toArray();
			//this.selectDisplacement(u_n_vector.toArray());
			this.assembleInternalForce(u_n_vector, r_int);

			IMatrix Kt = new Array2DMatrix(r_int.getN(), r_int.getN());
			this.assembleTangentStiffnessMatrix(Kt, u_n_vector);
			Vector rhs = r_ext.add(-1, r_int);
			Vector delta_u = this.solveLinearSystemEquation(Kt, rhs);
			convergence = delta_u.Norm2();
			System.out.println("convergence: " + convergence);
			u_n_vector = delta_u.add(1, u_n_vector);
		}
		;
		this.selectDisplacement(u_n_vector.toArray());
	}
}