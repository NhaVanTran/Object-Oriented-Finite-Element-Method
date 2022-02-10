package Models;

import FEM.*;
import inf.v3d.view.Viewer;
import java.lang.Math;

/**
 * The large Bottrop-Tetraeder.
 * 
 * M. Baitsch
 */
public class LargeTetraeder {

	public static Structure createStructure() {
		Structure struct = new Structure();
		double r1 = 457.2 / 2000;
		double r2 = 558.8 / 2000;
		double t = 10.0 / 1000;
		double a1 = Math.PI * (Math.pow(r1, 2) - Math.pow(r1 - t, 2));
		double a2 = Math.PI * (Math.pow(r2, 2) - Math.pow(r2 - t, 2));
		double e = 2.1e8;
		double l = 15.0;
		double u = Math.sqrt(3);
		double s = Math.sqrt(6);
		Constraint c1 = new Constraint(false, false, false);

		Node n1 = struct.addNode(0, l * u * 2 / 3, l * s * 4 / 3);
		Node n2 = struct.addNode(0, u * l, s * l);
		Node n3 = struct.addNode(-0.5 * l, u * l / 2, s * l);
		Node n4 = struct.addNode(0.5 * l, u * l / 2, s * l);
		Node n5 = struct.addNode(0, l * u * 4 / 3, l * s * 2 / 3);
		Node n6 = struct.addNode(-1 * l, u * l / 3, l * s * 2 / 3);
		Node n7 = struct.addNode(1 * l, u * l / 3, l * s * 2 / 3);
		Node n8 = struct.addNode(-0.5 * l, l * u * 5 / 6, l * s * 2 / 3);
		Node n9 = struct.addNode(0, u * l / 3, l * s * 2 / 3);
		Node n10 = struct.addNode(0.5 * l, l * u * 5 / 6, l * s * 2 / 3);
		Node n11 = struct.addNode(0, 2 * u * l, 0);
		Node n12 = struct.addNode(-2 * l, 0, 0);
		Node n13 = struct.addNode(2 * l, 0, 0);
		Node n14 = struct.addNode(-1 * l, u * l, 0);
		Node n15 = struct.addNode(1 * l, u * l, 0);
		Node n16 = struct.addNode(0, 0, 0);
		Node n17 = struct.addNode(-1.5 * l, u * 0.5 * l, 0);
		Node n18 = struct.addNode(1.5 * l, u * 0.5 * l, 0);
		Node n19 = struct.addNode(1 * l, 0, 0);
		Node n20 = struct.addNode(-1 * l, 0, 0);
		Node n21 = struct.addNode(1.5 * l, u * l / 6, s * l / 3);
		Node n22 = struct.addNode(-1.5 * l, u * l / 6, s * l / 3);
		Node n23 = struct.addNode(1 * l, l * u * 2 / 3, s * l / 3);
		Node n24 = struct.addNode(0.5 * l, u * l / 6, s * l / 3);
		Node n25 = struct.addNode(-0.5 * l, u * l / 6, s * l / 3);
		Node n26 = struct.addNode(-1 * l, l * u * 2 / 3, l * s / 3);
		Node n27 = struct.addNode(0.5 * l, u * 0.5 * l, 0);
		Node n28 = struct.addNode(-0.5 * l, u * 0.5 * l, 0);

		Force f2 = new Force(-45.4, -71.93, -90.09);
		Force f3 = new Force(139.42, 48.92, -205.98);
		Force f4 = new Force(-19.91, 19.36, -21.37);
		Force f5 = new Force(12.26, -125.010, -77.69);
		Force f6 = new Force(-49.52, 3.09, -88.08);
		Force f7 = new Force(13.15, -2.66, -6.7);
		Force f8 = new Force(134.01, 44.59, -109.7);
		Force f9 = new Force(-74.59, -25.81, -82.2);
		Force f10 = new Force(-146.26, 80.42, -144.34);
		Force f17 = new Force(-4.09, 74.42, -35.45);
		Force f20 = new Force(9.4, -51.35, -83.21);
		Force f22 = new Force(75.77, 27.4, -80.54);
		Force f25 = new Force(-63.55, 13.68, -45.79);
		Force f26 = new Force(-17.52, -64.18, -43.63);
		Force f28 = new Force(-57.27, -47.65, -106.9);

		n2.setForce(f2);
		n3.setForce(f3);
		n4.setForce(f4);
		n5.setForce(f5);
		n6.setForce(f6);
		n7.setForce(f7);
		n8.setForce(f8);
		n9.setForce(f9);
		n10.setForce(f10);
		n17.setForce(f17);
		n20.setForce(f20);
		n22.setForce(f22);
		n25.setForce(f25);
		n26.setForce(f26);
		n28.setForce(f28);

		n14.setConstraint(c1);
		n15.setConstraint(c1);
		n19.setConstraint(c1);
		n20.setConstraint(c1);

		struct.addElement(e, a1, 0, 1);
		struct.addElement(e, a1, 0, 3);
		struct.addElement(e, a1, 0, 2);
		struct.addElement(e, a1, 1, 2);
		struct.addElement(e, a1, 1, 3);
		struct.addElement(e, a1, 2, 3);
		struct.addElement(e, a1, 2, 7);
		struct.addElement(e, a1, 2, 8);
		struct.addElement(e, a1, 2, 5);
		struct.addElement(e, a1, 5, 7);
		struct.addElement(e, a1, 7, 8);
		struct.addElement(e, a1, 5, 8);
		struct.addElement(e, a1, 1, 4);
		struct.addElement(e, a1, 1, 9);
		struct.addElement(e, a1, 1, 7);
		struct.addElement(e, a1, 7, 4);
		struct.addElement(e, a1, 4, 9);
		struct.addElement(e, a1, 7, 9);
		struct.addElement(e, a1, 3, 9);
		struct.addElement(e, a1, 3, 6);
		struct.addElement(e, a1, 3, 8);
		struct.addElement(e, a1, 8, 9);
		struct.addElement(e, a1, 9, 6);
		struct.addElement(e, a1, 8, 6);
		struct.addElement(e, a2, 4, 10);
		struct.addElement(e, a2, 4, 14);
		struct.addElement(e, a2, 4, 13);
		struct.addElement(e, a2, 13, 10);
		struct.addElement(e, a2, 10, 14);
		struct.addElement(e, a2, 14, 13);
		struct.addElement(e, a1, 5, 25);
		struct.addElement(e, a1, 5, 24);
		struct.addElement(e, a1, 5, 21);
		struct.addElement(e, a1, 21, 25);
		struct.addElement(e, a1, 25, 24);
		struct.addElement(e, a1, 24, 21);
		struct.addElement(e, a1, 6, 22);
		struct.addElement(e, a1, 6, 20);
		struct.addElement(e, a1, 6, 23);
		struct.addElement(e, a1, 23, 22);
		struct.addElement(e, a1, 22, 20);
		struct.addElement(e, a1, 20, 23);
		struct.addElement(e, a1, 25, 13);
		struct.addElement(e, a1, 25, 27);
		struct.addElement(e, a1, 25, 16);
		struct.addElement(e, a1, 16, 13);
		struct.addElement(e, a1, 13, 27);
		struct.addElement(e, a1, 16, 27);
		struct.addElement(e, a1, 22, 14);
		struct.addElement(e, a1, 22, 17);
		struct.addElement(e, a1, 22, 26);
		struct.addElement(e, a1, 26, 14);
		struct.addElement(e, a1, 14, 17);
		struct.addElement(e, a1, 17, 26);
		struct.addElement(e, a1, 21, 16);
		struct.addElement(e, a1, 21, 19);
		struct.addElement(e, a1, 21, 11);
		struct.addElement(e, a1, 11, 16);
		struct.addElement(e, a1, 16, 19);
		struct.addElement(e, a1, 11, 19);
		struct.addElement(e, a1, 24, 27);
		struct.addElement(e, a1, 24, 15);
		struct.addElement(e, a1, 24, 19);
		struct.addElement(e, a1, 19, 27);
		struct.addElement(e, a1, 27, 15);
		struct.addElement(e, a1, 15, 19);
		struct.addElement(e, a1, 23, 26);
		struct.addElement(e, a1, 23, 18);
		struct.addElement(e, a1, 23, 15);
		struct.addElement(e, a1, 15, 26);
		struct.addElement(e, a1, 26, 18);
		struct.addElement(e, a1, 18, 15);
		struct.addElement(e, a1, 20, 17);
		struct.addElement(e, a1, 20, 12);
		struct.addElement(e, a1, 20, 18);
		struct.addElement(e, a1, 18, 17);
		struct.addElement(e, a1, 17, 12);
		struct.addElement(e, a1, 18, 12);

		return struct;
	}

	public static void main(String[] args) {
		Viewer viewer = new Viewer();
		Structure struct = createStructure();
		Visualizer viz = new Visualizer(struct, viewer);

		viz.radiusScale(5);
		viz.constrainScale(3);
		viz.forceScale(4 * 3e-2);
		viz.forceRadiusScale(.15);
		viz.displacementScale(1e3);
		viz.elementForceScale(1e-2);

		viz.drawElement();
		viz.drawConstraints();
		viz.drawForce();

		viewer.setVisible(true);
		// struct.solveLinear();
		struct.solveArclengthControl(5);
		// struct.solveNewtonRaphson(1, 1);
		viz.drawDisplacement();// Draw Displacement after solving problem
		viz.drawElementForce();
		// struct.printResult();

	}
}
