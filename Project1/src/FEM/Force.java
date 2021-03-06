package FEM;

import inf.text.ArrayFormat;

public class Force {
	private double[] components = new double[3];

	// Attribute
	public Force(double r1, double r2, double r3) {
		this.components[0] = r1;
		this.components[1] = r2;
		this.components[2] = r3;
	}

	// Get component of force vector
	public double getComponent(int c) {
		return this.components[c];
	}

	// printing
	public void print() {
		System.out.println(ArrayFormat.format(components));
	}
}
