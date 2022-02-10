package FEM;

import inf.jlinalg.Vector3D;
import inf.text.ArrayFormat;

public class Node {
	// Attribute of nodes
	private int[] dofNumbers = new int[3];
	private Force force = null;
	private Constraint constraint = null;
	private Vector3D displacement = null;
	private Vector3D position = null;

	/* Position of Node */
	public Node(double x1, double x2, double x3) {
		this.position = new Vector3D(x1, x2, x3);
	}

	public void setConstraint(Constraint c) {
		this.constraint = c;
	}

	public Constraint getConstraint() {
		return this.constraint;
	}

	public void setForce(Force f) {
		this.force = f;
	}

	public Force getForce() {
		return this.force;
	}

	public int enumerateDOFs(int start) {
		if (this.getConstraint() == null) { // if Node are not to be set
											// constraint, count DOF
			for (int i = 0; i < 3; i++) {
				this.dofNumbers[i] = start++;
			}
		} else {
			for (int i = 0; i < 3; i++) {
				if (this.getConstraint().isFree(i) == true) {
					this.dofNumbers[i] = start++;// free => count DOF
				} else {
					this.dofNumbers[i] = -1; // fixed
				}
			}
		}
		return start;
	}

	public int[] getDOFNumbers() {
		return this.dofNumbers;
	}

	public Vector3D getPosition() {
		return this.position;
	}

	public void setDisplacement(double[] u) {
		this.displacement = new Vector3D(u[0], u[1], u[2]);
	}

	public Vector3D getDisplacement() {
		return this.displacement;
	}

	public void print() {
		System.out.println(position.getX1() + "   " + position.getX2() + "   "
				+ position.getX3());
	}
}
