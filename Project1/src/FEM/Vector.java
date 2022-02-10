package FEM;

public class Vector {

	private double[] elements;

	public Vector(int n) {
		this.elements = new double[n];
	}

	public Vector(double... x) {
		this.elements = x;
	}

	public void assign(Vector a) {
		for (int i = 0; i < a.getN(); i++) {
			this.elements[i] = a.get(i);
		}
	}

	public double[] toArray() {
		double[] a = new double[this.getN()];
		for (int i = 0; i < a.length; i++) {
			a[i] = this.elements[i];
		}
		return a;
	}

	public int getN() {
		return this.elements.length;
	}

	public void fill(double v) {
		for (int i = 0; i < this.getN(); i++) {
			this.elements[i] = v;
		}
	}

	public void set(int idx, double v) {
		this.elements[idx] = v;
	}

	public double get(int idx) {
		return this.elements[idx];
	}

	public void print(String l) {

		if (getN() == 0) {
			System.out.print(l + " [];\n");
		} else if (getN() == 1) {
			System.out.print(l + " [ " + this.elements[0] + " ];\n");
		} else {
			System.out.print(l + " [ " + this.elements[0] + ", ");
			for (int i = 1; i < this.getN() - 1; i++) {
				System.out.print(this.elements[i] + ", ");
			}
			System.out.print(this.elements[this.getN() - 1] + " ]; \n");
		}

	}

	public Vector multiply(double alpha) {
		Vector z = new Vector(this.getN());
		for (int i = 0; i < this.getN(); i++) {
			z.set(i, this.elements[i] * alpha);
		}
		return z;
	}

	public Vector add(double alpha, Vector v) {
		Vector z = new Vector(this.getN());
		for (int i = 0; i < this.getN(); i++) {
			z.set(i, this.elements[i] + alpha * (v.get(i)));
		}
		return z;
	}

	public double scalarProduct(Vector v) {
		double s = 0;
		for (int i = 0; i < getN(); i++) {
			s += this.elements[i] * v.elements[i];
		}
		return s;
	}

	public double Norm2() {
		return Math.sqrt(this.scalarProduct(this));
	}

}
