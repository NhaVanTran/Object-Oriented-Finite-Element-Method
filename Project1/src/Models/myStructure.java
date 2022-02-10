package Models;

import inf.v3d.view.Viewer;

import java.util.ArrayList;

import FEM.Constraint;
import FEM.Force;
import FEM.Node;
import FEM.Structure;
import FEM.Visualizer;

public class myStructure {

	/**
	 * OOFEM @ Tran Van Nha CompEng 2012 M-08
	 */

	public static Structure createStructure() {
		Structure struct = new Structure();
		double r = 457.2 / 2000;
		double t = 10.0 / 1000;
		double a = Math.PI * (Math.pow(r, 2) - Math.pow(r - t, 2));
		double e = 2.1e11;
		Constraint c1 = new Constraint(false, false, false);
		Force f = new Force(1e-3, -10e3, 0);

		double Lz = 1000;// Length of roof in z direction
		double R1 = 300;// Smaller radius of roof
		double R2 = 400;// Bigger radius of roof
		int nz = 6;// element in z direction
		int nphi = 4;
		double dz = Lz / nz;
		double dphi = Math.PI / (2 * nphi);
		double[] phi = new double[nphi + 1];
		for (int i = 0; i < phi.length; i++) {
			phi[i] = Math.PI / 4 + dphi * i;
		}
		// create nodes
		for (int i = 0; i < nz + 1; i++) {
			for (int j = 0; j < nphi + 1; j++) {
				double x1 = R1 * Math.cos(phi[j]);
				double x2 = R1 * Math.sin(phi[j]);
				double x3 = i * dz;
				struct.addNode(x1, x2, x3);

				double X1 = R2 * Math.cos(phi[j]);
				double X2 = R2 * Math.sin(phi[j]);
				double X3 = i * dz;
				struct.addNode(X1, X2, X3);

			}
		}
		// struct.addNode(0,0, 0);
		// apply BCs
		for (int i = 0; i < nz + 1; i++) {
			struct.getNode(i * 2 * (nphi + 1)).setConstraint(c1);
			struct.getNode(i * 2 * (nphi + 1) + 2 * nphi).setConstraint(c1);
		}
		// Set force for odd node
		for (int i = 1; i < 2 * (nphi + 1) * (nz + 1); i = i + 2) {
			struct.getNode(i).setForce(f);
		}

		// create elements
		for (int i = 0; i < nphi; i++) {
			for (int j = 0; j < nz; j++) {
				int p = j * 2 * (nphi + 1) + 2 * i;
				struct.addElement(e, a, p, p + 2);
				struct.addElement(e, a, p, p + 3);
				struct.addElement(e, a, p, p + 1);
				struct.addElement(e, a, p + 1, p + 3);
				struct.addElement(e, a, p, p + 2 * (nphi + 1));

				struct.addElement(e, a, p, p + 2 * (nphi + 1) + 3);

				struct.addElement(e, a, p + 1, (p + 1) + 2 * (nphi + 1));

				struct.addElement(e, a, p + 2 * (nphi + 1), p + 2 * (nphi + 1)
						+ 2);
				struct.addElement(e, a, p + 2 * (nphi + 1), p + 2 * (nphi + 1)
						+ 3);
				struct.addElement(e, a, p + 2 * (nphi + 1), p + 2 * (nphi + 1)
						+ 1);
				struct.addElement(e, a, p + 2 * (nphi + 1) + 1, p + 2
						* (nphi + 1) + 3);

			}
		}
		for (int j = 0; j < nz; j++) {
			int p = j * 2 * (nphi + 1) + 2 * nphi;
			struct.addElement(e, a, p, p + 2 * (nphi + 1));
			struct.addElement(e, a, p + 1, (p + 1) + 2 * (nphi + 1));
			struct.addElement(e, a, p, p + 1);
			struct.addElement(e, a, p + 2 * (nphi + 1), p + 1 + 2 * (nphi + 1));
		}

		// return the new structure
		return struct;
	}

	public static void main(String[] args) {

		Viewer viewer = new Viewer();
		Structure struct = createStructure();
		Visualizer viz = new Visualizer(struct, viewer);

		viz.radiusScale(50);
		viz.radiusNodeScale(10);
		viz.constrainScale(30);
		viz.forceScale(1e-2);
		viz.forceRadiusScale(2);
		viz.displacementScale(1e4);
		viz.elementForceScale(3e-12);

		viz.drawNode();
		viz.drawElement();
		viz.drawConstraints();
		viz.drawForce();

		viewer.setVisible(true);

		struct.solveLinear();
		// struct.SolveNonLinearArcLengthControl(5);
		viz.drawDisplacement();// Draw Displacement after solving problem
		viz.drawElementForce();
		struct.printResult();

	}

}
