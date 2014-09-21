class Point extends PVector {
  public Point() {
    this(0, 0);
  }
  public Point(float x, float y) {
    super(x, y);
  }
}

class Rectangle {
  public float x, y, width, height;
  public Rectangle() {
    this(0, 0, 0, 0);
  }
  public Rectangle(float x, float y, float width, float height) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }
}

class Particle {
  public Point p0; // current
  public Point p1; // next
  public Point v;
  public float mass = 1.0;
  public Particle() {
    this(0, 0, 0, 0);
  }
  public Particle(float x, float y) {
    this(x, y, 0, 0);
  }
  public Particle(float x, float y, float vx, float vy) {
    p0 = new Point(x, y);
    p1 = new Point(x, y);
    v = new Point(vx, vy);
  }
  public void set(float x, float y, float vx, float vy) {
    setPosition(x, y);
    setNextPosition(x, y);    
    setVelocity(vx, vy);
  }
  public void setPosition(float x, float y) {
    p0.set(x, y);
  }
  public void setNextPosition(float x, float y) {
    p1.set(x, y);
  }
  public void setVelocity(float vx, float vy) {
    v.set(vx, vy);
  }
  public float x() {
    return p0.x;
  }
  public float y() {
    return p0.y;
  }
  public float x1() {
    return p1.x;
  }
  public float y1() {
    return p1.y;
  }
  public float vx() {
    return v.x;
  }
  public float vy() {
    return v.y;
  }
}

class Constraint {
  // cardinality
  public int n = 0;
  
  // indices
  // 1...n
  public int[] indices;
  
  // positions
  // 1...n
  public Particle[] positions;
  
  // stiffness
  // 0...1
  public float k = 1.0f;

  // equality or inequality
  //   equality: eval() == 0
  // inequality: eval() >= 0
  public boolean isEquality = true;
  
  public Constraint(int... indices) {
    this.n = indices.length;
    this.indices = indices;
    this.positions = new Particle[n];
  }
  
  // evaluation function
  // R^{3n} -> R
  public float eval(int iterations) {
    return 0;
  }
  
  public void setPosition(int i, Particle p) {
    positions[i] = p;
  }
  
  public Particle getPosition(int i) {
    return positions[i];
  }
  
  protected float _k(int ns) {
    return 1.0f - pow(1.0f - k, 1.0f / ns);
  }
}

class DistanceConstraint extends Constraint {
  public float distance; 
  public DistanceConstraint(int i0, int i1, float distance) {
    super(i0, i1);
    this.distance = distance;
  }  
  public float eval(int iterations) {   
    Point p1 = positions[0].p1;
    Point p2 = positions[1].p1;
//    Point p = new Point(p1.x - p2.x, p1.y - p2.y);
    float d = p1.dist(p2);
    float c = d - distance;
    
    if(!isEquality && c >= 0) {
      return c;
    }
    
    float m1 = positions[0].mass;
    float m2 = positions[1].mass;
    float w1 = m1 > 0 ? 1.0 / m1 : 0;
    float w2 = m2 > 0 ? 1.0 / m2 : 0;

    if(w1 == 0 && w2 == 0) {
      return c;
    }
    
    float s = c / (w1 + w2);
    PVector n = PVector.div(PVector.sub(p1, p2), d);
    
    Point dp1 = new Point(-s * w1 * n.x, -s * w1 * n.y);
    Point dp2 = new Point( s * w2 * n.x,  s * w2 * n.y);
    
    float dt = 1.0 / 30; // TODO
    float t = _k(iterations) * dt;
    p1.set(p1.x + t * dp1.x, p1.y + t * dp1.y);
    p2.set(p2.x + t * dp2.x, p2.y + t * dp2.y);
    
    return c;
  }
}

// C(p) = dot((p - q), n)
class CollisionConstraint extends Constraint {
  Point query;
  Point normal;
  public CollisionConstraint(int i, Point query, Point normal) {
    super(i);
    this.query = query;
    this.normal = normal;
  }
  public float eval(int iterations) {
    Point p = positions[0].p1;
    Point q = query;
    Point qp = new Point(p.x - q.x, p.y - q.y);
    
    // C(p) = dot((p-q), n)
    float c = PVector.dot(qp, normal);
    
    if(!isEquality && c >= 1.0f) {
      return c;
    }
    
    if(positions[0].mass <= 0) {
      return c;
    }
    
    // delta p = -s w nabla C(p)
    Point dp = new Point(-qp.x, -qp.y);
    
    float dt = 1.0f / 30; // TODO
    float t = _k(iterations) * dt;
    p.set(p.x + t * dp.x, p.y + t * dp.y);
    
    return c;
  }
}

class BendConstraint extends Constraint {
  class Result {
    public PVector n0, n1;
    public float d, phi;
  }
  public float phi0;
  private Result _res;
  PVector _cpp12, _cpp13, _cpn11, _cnp01, _cpn10, _cnp11, _cpn21, _cnp02, _cpn30, _cnp13;
  public BendConstraint(int i0, int i1, int i2, int i3, Particle pa0, Particle pa1, Particle pa2, Particle pa3) {
    super(i0, i1, i2, i3);
    _res = new Result();
    _cpp12 = new PVector(); 
    _cpp13 = new PVector(); 
    _cpn11 = new PVector(); 
    _cnp01 = new PVector(); 
    _cpn10 = new PVector(); 
    _cnp11 = new PVector(); 
    _cpn21 = new PVector(); 
    _cnp02 = new PVector(); 
    _cpn30 = new PVector(); 
    _cnp13 = new PVector();
    // calc phi0
    Point p0 = pa0.p0, p1 = pa1.p0, p2 = pa2.p0, p3 = pa3.p0;
    phi0 = _phi(p0, p1, p2, p3).phi;
    println("phi9=" + phi0);
  }
  public float eval(int interations) {
    Point p0 = positions[0].p1;
    Point p1 = positions[1].p1;
    Point p2 = positions[2].p1;
    Point p3 = positions[3].p1;
    
    Result res = _phi(p0, p1, p2, p3);
    float d = res.d;
    float c = res.phi - phi0;
    
    if(!isEquality && c >= 1.0f) {
      return c;
    }
    
    PVector n0 = _res.n0, n1 = _res.n1;
    
    PVector.cross(p1, p2, _cpp12);
    PVector.cross(p1, p3, _cpp13);
    PVector.cross(p1, n1, _cpn11);
    PVector.cross(n0, p1, _cnp01);
    PVector.cross(p1, n0, _cpn10);
    PVector.cross(n1, p1, _cnp11);
    PVector.cross(p2, n1, _cpn21);
    PVector.cross(n0, p2, _cnp02);
    PVector.cross(p3, n0, _cpn30);
    PVector.cross(n1, p3, _cnp13);
    
    float d12 = _cpp12.mag();
    float d13 = _cpp13.mag();
    
    PVector q2 = PVector.div(PVector.add(_cpn11, PVector.mult(_cnp01, d)), d12);
    PVector q3 = PVector.div(PVector.add(_cpn10, PVector.mult(_cnp11, d)), d13);
    PVector q20 = PVector.div(PVector.add(_cpn21, PVector.mult(_cnp02, d)), d12);
    PVector q21 = PVector.div(PVector.add(_cpn30, PVector.mult(_cnp13, d)), d13);
    PVector q1 = PVector.add(PVector.mult(q20, -1), PVector.mult(q21, -1));
    PVector q0 = PVector.sub(PVector.sub(PVector.mult(q1, -1), q2), q3);
    
    // TODO: calc delta p_i

    return c;
  }
  Result _phi(PVector p0, PVector p1, PVector p2, PVector p3) {
    PVector n0 = _res.n0 = _normal(p0, p1, p2);
    PVector n1 = _res.n1 = _normal(p0, p1, p3);
    float d = _res.d = PVector.dot(n0, n1);
    _res.phi = acos(d);
    return _res;
  }
  PVector _normal(PVector p0, PVector p1, PVector p2) {
    PVector c0 = new PVector();
    PVector.cross(PVector.sub(p1, p0), PVector.sub(p2, p0), c0); 
    float d0 = c0.mag();
    if(d0 == 0) {
      return new PVector(0, 0, 0);
    }
    return PVector.div(c0, d0);
  }
}

class PObject {
  
}

class PBox extends PObject {
  public Rectangle bounds;
  public PBox(float x, float y, float width, float height) {
    bounds = new Rectangle(x, y, width, height);
  }
}

float _lastUpdateTime;
boolean _isPaused = false;

ArrayList<Particle> _particles;
ArrayList<Point> _forces; // external force
ArrayList<Constraint> _constraints;
ArrayList<PObject> _objects;
int _solverIterations = 4;
int _updateIterations = 2;

void setup() {
  size(500, 500);
  frameRate(30);
  
  _lastUpdateTime = millis();
  
  _particles = new ArrayList<Particle>();
  _forces = new ArrayList<Point>();
  _constraints = new ArrayList<Constraint>();
  
//  Particle p = new Particle();
//  p.setPosition(width / 2, height / 2);
//  p.setVelocity(0, 13);
//  p.mass = 0;
//  _particles.add(p);
//  
//  Particle p2 = new Particle();
//  p2.setPosition(p.x() + 100, p.y());
//  p2.setVelocity(0, 13);
//  _particles.add(p2);
//  
//  DistanceConstraint dc = new DistanceConstraint(0, 1, 50);
//  dc.k = 0.5f;
//  _constraints.add(dc);
//  
//  _forces.add(new Point(10, 0));

  float w = 100, h = 100;
  int segX = 10, segY = 10;
  float tx = (width - w) / 2, ty = (height - h) / 2;
  for(int iy = 0; iy <= segY; ++ iy) {
    float y = ty + iy * h / segY;
    for(int ix = 0; ix <= segX; ++ ix) {
      float x = tx + ix * w / segX;
      Particle p = new Particle(x, y, 0, 0);
      p.mass = iy == 0 ? 0 : 1.0f;
      _particles.add(p);
    }
  }
  
  int step = segX + 1;
  float dw = w / segX, dh = h / segY, dr = sqrt(dw * dw + dh * dh);
  for(int iy = 0; iy < segY; ++ iy) {  
    for(int ix = 0; ix < segX; ++ ix) {
      int i = iy * step + ix;
      int i0 = i, i1 = i + 1, i2 = i + step, i3 = i2 + 1;
      _constraints.add(new DistanceConstraint(i3, i2, dw));
      _constraints.add(new DistanceConstraint(i2, i0, dh));
      _constraints.add(new DistanceConstraint(i0, i3, dr));
      _constraints.add(new DistanceConstraint(i3, i0, dr));
      _constraints.add(new DistanceConstraint(i0, i1, dw));
      _constraints.add(new DistanceConstraint(i1, i3, dh));      
    }
  }
 
//  _particles.add(new Particle(tx, ty, 0, 0));
//  _particles.add(new Particle(tx + w, ty, 0, 0));
//  _particles.add(new Particle(tx, ty + h, 0, 0));
//  _particles.add(new Particle(tx + w, ty + h, 0, 0));
//  
//  _particles.get(0).mass = _particles.get(1).mass = 0;
//   
//  float t = sqrt(w * w + h * h);
//  _constraints.add(new DistanceConstraint(3, 2, w));
//  _constraints.add(new DistanceConstraint(2, 0, h));
//  _constraints.add(new DistanceConstraint(0, 3, t));
//  _constraints.add(new DistanceConstraint(3, 0, t));
//  _constraints.add(new DistanceConstraint(0, 1, w));
//  _constraints.add(new DistanceConstraint(1, 3, h));
  
//  _constraints.add(new BendConstraint(3, 0, 2, 1,
//                                      _particles.get(3),
//                                      _particles.get(0),
//                                      _particles.get(2),
//                                      _particles.get(1)));
  
  for(int i = 0, len = _constraints.size(); i < len; ++ i) {
    Constraint c = _constraints.get(i);
    c.k = 1.0f;
    c.isEquality = true;
  }
  
  _forces.add(new Point(0, 100));
}

void draw() {
  // update
  update();

  // draw  
  background(0);
  
  for(int i = 0; i < _particles.size(); ++ i) {
    Particle p = _particles.get(i);
    drawParticle(p);
//    if(i == 0) {
//      text("p=(" + p.p0 + ") \nv=(" + p.v + ")", 10, 10);
//    }
  }
  
  noFill();
  stroke(50, 150, 255);
  for(int i = 0, len = _constraints.size(); i < len; ++ i) {
    Constraint c = _constraints.get(i);
    if(c instanceof DistanceConstraint) {
      DistanceConstraint dc = (DistanceConstraint) c;
      Particle p0 = _particles.get(dc.indices[0]);
      Particle p1 = _particles.get(dc.indices[1]);
      line(p0.x(), p0.y(), p1.x(), p1.y());
    }
  }
  
  fill(255);
  stroke(255);
  text("" + frameRate, 10, height - 10);
}

void update() {
//  if(_isPaused) {
//    _lastUpdateTime = millis();
//    return;
//  }
//  _isPaused = true;
 
  float time = millis();

  float dt = (time - _lastUpdateTime) / 1000.0f;
  _lastUpdateTime = time;
  
  dt = dt / (float) _updateIterations;
  for(int j = 0; j < _updateIterations; ++ j) { 
  
  Point force = new Point(0, 0);
  for(int i = 0; i < _forces.size(); ++ i) {
    Point f = _forces.get(i);
    force.x += f.x;
    force.y += f.y;
  }

  updateParticles(_particles, force, dt);
  }
}

void updateParticles(ArrayList<Particle> particles, Point force, float dt) {
  // forall vertices i do vi ← vi + ∆t wi fext (xi)
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    float mass = p.mass;
    if(mass > 0) {
      float vx = p.vx() + dt * force.x / mass;
      float vy = p.vy() + dt * force.y / mass;
      p.setVelocity(vx, vy);
    }
  }
  // dampVelocities(v1 , . . . , vN )

  // forall vertices i do pi ← xi + ∆t vi
  Point[] ps = new Point[particles.size()];
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    if(p.mass > 0) {
      float x = p.x() + dt * p.vx();
      float y = p.y() + dt * p.vy();
      p.setNextPosition(x, y);
    } else {
      p.setNextPosition(p.x(), p.y());
    }
  }
  // forall vertices i do generateCollisionConstraints(xi → pi)
  generateCollisionConstraints(particles, ps);
  
  // loop solverIterations times
  //   projectConstraints(C1,...,CM+Mcoll ,p1,...,pN)
  // endloop
  for(int i = 0; i < _solverIterations; ++ i) {
    projectConstraints(_constraints, _particles);
  }
  
  // forall vertices i
  //   vi ←(pi−xi)/∆t
  //   xi ← pi
  // endfor
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    float x = p.p1.x;
    float y = p.p1.y;
    float vx = (x - p.x()) / dt;
    float vy = (y - p.y()) / dt;
    p.set(x, y, vx, vy);
  }
  
  // velocityUpdate(v1 , . . . , vN )
  
}

ArrayList<Constraint> generateCollisionConstraints(ArrayList<Particle> xs, Point[] ps) {
  return null;
}

void projectConstraints(ArrayList<Constraint> constraints, ArrayList<Particle> particles) {
  for(int i = 0; i < constraints.size(); ++ i) {
    Constraint c = constraints.get(i);
    for(int j = 0; j < c.n; ++ j) {
      c.setPosition(j, particles.get(c.indices[j]));
    }
    float r = c.eval(_solverIterations);
//    for(int j = 0; j < c.n; ++ j) { 
//      positions[c.indices[j]].set(c.getPosition(j));      
//    }
  }
}

//void updateParticle(Particle p, Point f, float dt) { 
//  // prediction
//  float x = p.x() + dt * p.vx();
//  float y = p.y() + dt * p.vy();
//  // position correction
//  if(y > 400) { y = 400; }
//  if(x > 400) { x = 400; }
//  // velocity update
//  float vx = (x - p.x()) / dt;
//  float vy = (y - p.y()) / dt;  
//  // velocity correction
//  vx = vx + dt * f.x / p.mass;
//  vy = vy + dt * f.y / p.mass;
//  // set property
//  p.set(x, y, vx, vy);
//}

void mousePressed() {
  _isPaused = false;
}

void mouseDragged() {
//  int i = _particles.size() - 1;
//  _particles.get(i).setPosition(mouseX, mouseY);
  float dx = mouseX - width / 2, dy = mouseY - height / 2;
  _forces.get(0).set(dx * 0.1, dy * 0.1);
}

void drawParticle(Particle p) {
  float x = p.x(), y = p.y();
  fill(255);
  noStroke();
  ellipse(x, y, 4, 4);
  noFill();
  stroke(255,128,0);
  line(x, y, x + p.vx(), y + p.vy());
}
