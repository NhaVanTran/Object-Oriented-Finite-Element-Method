package FEM;

import inf.text.ArrayFormat;

public class Constraint {

	private boolean[] free = new boolean[3];

	/* Constraint */
	public Constraint(boolean b1, boolean b2, boolean b3) {
		this.free[0] = b1;	// In X direction
		this.free[1] = b2;	// In Y direction
		this.free[2] = b3;	// In Z direction
	}

	/* Checking whether boundary Condition is free or not */
	public boolean isFree(int c) {
		return this.free[c];
	}

	// Printing
	public void print() {
		System.out.println(ArrayFormat.format(free));
	}
}
