package Test;

import inf.v3d.obj.Cone;
import inf.v3d.obj.Sphere;
import inf.v3d.view.Viewer;

public class drawCone {

	/**
	 * Test how to set direction for a cone
	 */
	public static void main(String[] args) {

		Viewer v = new Viewer();
		double a[] = { 11, 5, 3 };
		double b[] = { a[0] + 10, a[1], a[2] };
		double direc[] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };// Direction
																	// from a to
																	// b
		double lcone = Math.sqrt(Math.pow((a[0] - b[0]), 2)
				+ Math.pow((a[1] - b[1]), 2) + Math.pow((a[2] - b[2]), 2));

		Cone cone = new Cone();
		cone.setCenter(a[0], a[1], a[2]);
		cone.setHeight(lcone);
		cone.setRadius(lcone / 2);
		cone.setDirection(direc[0], direc[1], direc[2]);
		cone.setColor("green");

		Cone cone2 = new Cone();
		cone2.setCenter(a[0], a[1], a[2]);
		cone2.setHeight(lcone);
		cone2.setRadius(lcone / 2);
		cone2.setDirection(-direc[0], -direc[1], -direc[2]);
		cone2.setColor("blue");

		// Draw node a and node b to understand how to set direction of cone
		// In direction from a to b : d= b-a
		// otherwise : d= a-b
		Sphere sa = new Sphere(a[0], a[1], a[2]);
		Sphere sb = new Sphere(b[0], b[1], b[2]);
		sa.setRadius(lcone / 10);
		sb.setRadius(lcone / 10);
		sa.setColor("red");
		sb.setColor("cyan");
		v.addObject3D(sa);
		v.addObject3D(sb);
		v.addObject3D(cone);
		 v.addObject3D(cone2);

		v.setVisible(true);

		
	}

}
