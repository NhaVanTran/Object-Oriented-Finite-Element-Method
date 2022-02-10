package Test;

import inf.jlinalg.Vector3D;
import inf.v3d.obj.Cone;
import inf.v3d.obj.CylinderSet;
import inf.v3d.obj.Sphere;
import inf.v3d.view.Viewer;
import inf.v3d.obj.Cylinder;

public class TestVisualizer {

	public static void main(String[] args) {
		Viewer v = new Viewer();
		CylinderSet cs = new CylinderSet();
		Vector3D p1 = new Vector3D(1, 0, 0);
		Vector3D p2 = new Vector3D(10, 0, 0);
		double r = 1;
		cs.addCylinder(p1.toArray(), p2.toArray(), r);
		v.addObject3D(cs);

		// Draw constraint

		double a[] = { 11, 0, 0 };

		Cone cone = new Cone();
		cone.setCenter(a[0], a[1], a[2]);
		cone.setHeight(1);
		cone.setRadius(2);
		cone.setDirection(-1, 0, 0);
		v.addObject3D(cone);

		v.setVisible(true);

		// double a[]={1,2,3};
		// double b[]={4,5,6};
		// double c=a-b;
		// System.out.println(c);
	}

}
