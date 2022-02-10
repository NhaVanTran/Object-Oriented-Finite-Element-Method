package FEM;

import java.util.ArrayList;

import inf.jlinalg.Vector3D;
import inf.v3d.obj.Arrow;
import inf.v3d.obj.Cone;
import inf.v3d.obj.CylinderSet;
import inf.v3d.obj.PolygonSet;
import inf.v3d.obj.Sphere;
import inf.v3d.view.Viewer;

public class Visualizer {

	private Structure struct = new Structure();
	private Viewer viewer = new Viewer();
	private double displacementScale;
	private double forceScale;
	private double forceRadiusScale;
	private double radiusScale;
	private double radiusNodeScale;
	private double constraintScale;
	private double elementForceScale;

	public Visualizer(Structure struct, Viewer v) {
		this.struct = struct;
		this.viewer = v;
	}

	public void displacementScale(double s) {
		this.displacementScale = s;
	}

	public void forceScale(double s) {
		this.forceScale = s;
	}

	public void forceRadiusScale(double s) {
		this.forceRadiusScale = s;
	}

	public void radiusScale(double s) {
		this.radiusScale = s;
	}

	public void radiusNodeScale(double s) {
		this.radiusNodeScale = s;
	}

	public void constrainScale(double s) {
		this.constraintScale = s;
	}

	public void elementForceScale(double s) {
		this.elementForceScale = s;
	}

	public void drawElement() {
		CylinderSet cs = new CylinderSet();
		// Number of element
		int nele = this.struct.getNumberOfElement();
		// draw Cylinder
		for (int i = 0; i < nele; i++) {
			double[] p1 = this.struct.getElement(i).getNode1().getPosition()
					.toArray();// Position of node 1 of element
			double[] p2 = this.struct.getElement(i).getNode2().getPosition()
					.toArray();// Position of node 2 of element
			double r = Math.sqrt(this.struct.getElement(i).getArea() / Math.PI);
			// Scale radius
			r = r * this.radiusScale;
			cs.addCylinder(p1, p2, r);
			// cs.setColor("gray");
			this.viewer.addObject3D(cs);
		}
	}

	public void drawNode() {
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
			//double[] p = this.struct.getNode(i).getPosition().toArray();
			Sphere sp = new Sphere( 
					this.struct.getNode(i).getPosition().getX1(),  
					this.struct.getNode(i).getPosition().getX2(), 
					this.struct.getNode(i).getPosition().getX3());
			double r = this.radiusNodeScale;
			sp.setRadius(r);
			this.viewer.addObject3D(sp);
		}

	}

	public void drawConstraints() {
		int nnode = struct.getNumberOfNodes(); // Number of nodes
		for (int i = 0; i < nnode; i++) {
			Node node = this.struct.getNode(i);
			double[] pnode = node.getPosition().toArray();// Position of node(i)
			if (node.getConstraint() != null) {
				for (int j = 0; j < 3; j++) {
					if (node.getConstraint().isFree(j) == false) {
						// Length of cone
						double lcone = this.constraintScale;
						double[] center = this.struct.getNode(i).getPosition()
								.toArray();
						// Center of cone that we want to draw is in the lower
						// position of node
						center[j] = center[j] - lcone;
						/*
						 * NOTE: In java "center = pnode" means center point to
						 * object that pnode pointed to. Then the
						 * "center[j] = center[j] - lcone;" will change the
						 * value of object. More details, pnode[j] will also be
						 * changed. There for the lower line will lead to a
						 * result that direction ={0,0,0}; => Cannot draw cone
						 */

						// Direction of cone from center to node
						double[] direction = { pnode[0] - center[0],
								pnode[1] - center[1], pnode[2] - center[2] };
						// Radius of cone
						double r = lcone / 3;
						// Draw cone
						Cone cone = new Cone();
						cone.setCenter(center[0], center[1], center[2]);
						cone.setHeight(lcone);
						cone.setRadius(r);
						cone.setDirection(direction[0], direction[1],
								direction[2]);
						cone.setColor("green");
						this.viewer.addObject3D(cone);
					}
				}
			}
		}
	}

	public void drawForce() {
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
			Node node = this.struct.getNode(i);
			if (node.getForce() != null) {
				double[] pnode = node.getPosition().toArray();
				// Point 2 of scaled Arrow
				double[] p2 = new double[3];
				double p2x = node.getForce().getComponent(0) * this.forceScale;
				double p2y = node.getForce().getComponent(1) * this.forceScale;
				double p2z = node.getForce().getComponent(2) * this.forceScale;
				p2[0] = pnode[0] + p2x;
				p2[1] = pnode[1] + p2y;
				p2[2] = pnode[2] + p2z;
				// Point 1 of scaled Arrow is the position of node
				Arrow a = new Arrow(pnode, p2);
				a.setColor("red");
				a.setRadius(this.forceRadiusScale);
				this.viewer.addObject3D(a);
			}
		}
	}

	public void drawDisplacement() {
		// Draw deformed structure
		CylinderSet cs = new CylinderSet();
		// draw Cylinder
		for (int i = 0; i < this.struct.getNumberOfElement(); i++) {
			// The idea to draw displacement is adding displacement which is
			// scaled to position of node. Then draw element after deformation
			Vector3D n1 = this.struct.getElement(i).getNode1().getPosition();
			Vector3D n2 = this.struct.getElement(i).getNode2().getPosition();
			Vector3D u1 = this.struct.getElement(i).getNode1()
					.getDisplacement();
			Vector3D u2 = this.struct.getElement(i).getNode2()
					.getDisplacement();
			u1 = u1.multiply(this.displacementScale);
			u2 = u2.multiply(this.displacementScale);

			double[] p1 = n1.add(u1).toArray();// Position of node 1 of element
			double[] p2 = n2.add(u2).toArray();// Position of node 2 of element
			double r = Math.sqrt(this.struct.getElement(i).getArea() / Math.PI);
			// Scale radius
			r = r * this.radiusScale;
			cs.addCylinder(p1, p2, r);
			cs.setColor("orange");
			this.viewer.addObject3D(cs);
		}
	}

	public void drawElementForce() {
		for (int i = 0; i < this.struct.getNumberOfElement(); i++) {
			// Firstly, we find the vector that is perpendicular to E1-uniaxial
			// element vector
			Element e = this.struct.getElement(i);
			Vector3D x1 = this.struct.getElement(i).getNode1().getPosition();
			Vector3D x2 = this.struct.getElement(i).getNode2().getPosition();
			Vector3D n1 = new Vector3D(1, 0, 0);
			Vector3D n2 = new Vector3D(0, 1, 0);

			Vector3D p = n1.vectorProduct(this.struct.getElement(i).getE1());// p=n1xE1
			if (p.norm2() == 0) {
				// In case p parallel to n1
				p = n2.vectorProduct(this.struct.getElement(i).getE1());// p=n1xE1
			}
			// s1=x1+elementForceScale*ElementForce*p;
			Vector3D s1 = x1.add(p.multiply(e.computeForce()).multiply(
					this.elementForceScale));
			Vector3D s2 = x2.add(p.multiply(e.computeForce()).multiply(
					this.elementForceScale));
			// Add 4 vertex to a polygon

			PolygonSet ps = new PolygonSet();

			ps.insertVertex(x1.getX1(), x1.getX2(), x1.getX3(),
					e.computeForce());
			ps.insertVertex(x2.getX1(), x2.getX2(), x2.getX3(),
					e.computeForce());
			ps.insertVertex(s2.getX1(), s2.getX2(), s2.getX3(),
					e.computeForce());
			ps.insertVertex(s1.getX1(), s1.getX2(), s1.getX3(),
					e.computeForce());
			ps.polygonComplete();

			ps.setColoringByData(true);
			this.viewer.addObject3D(ps);

		}
	}
}
