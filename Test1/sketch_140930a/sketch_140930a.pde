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
  public Point q; // query
  public Point n; // normal
  public boolean isStatic;
  public CollisionConstraint(int i, Point q, Point n, boolean isStatic) {
    super(i);
    this.q = new Point(q.x, q.y);
    this.n = new Point(n.x, n.y);
    this.isStatic = isStatic;
    this.isEquality = false;
  }
  public float eval(int iterations) {
    Point p = positions[0].p1;
    Point qp = new Point(p.x - q.x, p.y - q.y);
    
    // C(p) = dot((p-q), n)
    float c = PVector.dot(qp, n);
    
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
  public float eval(int iterations) {
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

// C(p_1...p_n) = (sigma (p_i x p_i+1) / 2) - k_p A0
// C(p_1...p_n) = A - k_p A0
class Pressure2DConstraint extends Constraint {
  public float A0;
  public float kp; // k_pressure
  public Pressure2DConstraint(int[] indices, Particle[] points) {
    super(indices);
    A0 = _area(points);
    kp = 1.0f;
    println("A0=" + A0);
  }
  public float eval(int iterations) {
    
    float A = _area(positions);
    float c = A - kp * A0;
    
    if(!isEquality && c >= 1.0f) {
      return c;
    }
    
    int size = positions.length;
    if(size <= 2) {
      return c;
    }
    
    PVector[] npc = new PVector[size]; // nablas
    float[] ws = new float[size];
    float denom = 0;
    Point p0 = positions[size-1].p1;
    for(int i = 0; i < size; ++ i) {
      Point p1 = positions[i].p1;
      Point p2 = positions[(i+1)%size].p1;
      npc[i] = new PVector(p2.y - p0.y, -p2.x + p0.x);
      p0 = p1;
      //
      float mass = positions[i].mass;
      ws[i] = mass <= 0 ? 0 : 1.0f / mass;
      denom += ws[i] * npc[i].magSq();
    }
    
    float s = c / denom;
    
    float dt = 1.0f / 30; // TODO
    float t = _k(iterations) * dt;

    // delta p
    for(int i = 0; i < size; ++ i) {
      float dx = -s * ws[i] * npc[i].x;
      float dy = -s * ws[i] * npc[i].y;
      Point p = positions[i].p1;
      p.set(p.x + t * dx, p.y + t * dy);
    }
    
    return c;
  }
  float _area(Particle[] points) {
    int size = points.length;
    if(size <= 2) {
      return 0;
    }
    
    float area = 0;
    
    Point p0 = points[size-1].p1;
    for(int i = 0; i < size; ++ i) {
      Point p1 = points[i].p1;
      area += (p0.x * p1.y - p0.y * p1.x);
      p0 = p1;
    }
    
    return abs(area / 2);
  }
}

class Bend2DConstraint extends Constraint {
  class Result {
    PVector n0, n1;
    float d, phi;
  }
  public float phi0;
  private Result _res;
  public Bend2DConstraint(int i0, int i1, int i2, Particle pa0, Particle pa1, Particle pa2) {
    super(i0, i1, i2);
    _res = new Result();
    phi0 = _phi(pa0.p1, pa1.p1, pa2.p1).phi;
    println("phi0=" + (phi0 * 180 / PI));
  }
  Result _phi(PVector p0, PVector p1, PVector p2) {
    _res.n0 = PVector.sub(p1, p0);
    _res.n0.normalize();
    _res.n1 = PVector.sub(p2, p1);
    _res.n1.normalize();
    _res.d = PVector.dot(_res.n0, _res.n1);
    _res.phi = acos(_res.d);
    return _res;
  }  
}

// C(p) = |p-p0|
class SpringConstraint extends Constraint {
  public Point p0;
  public SpringConstraint(int i, Point p0) {
    super(i);
    this.p0 = new Point(p0.x, p0.y);
  }
  public float eval(int iterations) {
    Point p1 = positions[0].p1;
    float d = p1.dist(p0);
    float c = d;
    
    if(!isEquality && c >= 0) {
      return c;
    }
    
    float m1 = positions[0].mass;
    float w1 = m1 > 0 ? 1.0 / m1 : 0;

    if(w1 == 0) {
      return c;
    }
    
    PVector n = PVector.div(PVector.sub(p1, p0), d);
    
    Point dp1 = new Point(-(p1.x - p0.x), -(p1.y - p0.y));
    
    float dt = 1.0 / 30; // TODO
    float t = _k(iterations) * dt;
    p1.set(p1.x + t * dp1.x, p1.y + t * dp1.y);
    
    return c;
  }
}

class RigidBody { 
  public boolean inside(Point p) {
    return false;
  }
  public Point normal(Point p) {
    return null;
  } 
  public Point closestPoint(Point p) {
    return null;
  }
  public Point intersection(Point p0, Point p1) {
    return null;
  }
}

class BoxBody extends RigidBody {
  public Rectangle bounds;
  public BoxBody(float x, float y, float width, float height) {
    bounds = new Rectangle(x, y, width, height);
  } 
}

class CircleBody extends RigidBody {
  public Point center;
  public float radius;
  public CircleBody(float x, float y, float radius) {
    this.center = new Point(x, y);
    this.radius = radius;
  }
  public boolean inside(Point p) {
    float dx = p.x - center.x;
    float dy = p.y - center.y;
    return dx * dx + dy * dy <= radius * radius;
  }
  public Point normal(Point p) {
    float dx = p.x - center.x;
    float dy = p.y - center.y;
    float d = sqrt(dx * dx + dy * dy);
    return d == 0 ? new Point(1, 0) : new Point(dx / d, dy / d);
  }
  public Point closestPoint(Point p) {
    Point n = normal(p);    
    return new Point(center.x + n.x * radius, center.y + n.y * radius);
  }
  public Point intersection(Point p0, Point p1) {
    float vx0 = p1.x - p0.x, vy0 = p1.y - p0.y;    
    if(vx0 == 0 && vy0 == 0) {
      return null;
    }
    
//    println("v", vx0, vy0);
    
    float vx1 = center.x - p0.x, vy1 = center.y - p0.y;
    float d0 = sqrt(vx0 * vx0 + vy0 * vy0);
    float cross = (vx0 * vy1 - vy0 * vx1) / d0;
    float h = abs(cross);
    
    if(h > radius) {
      return null;
    }
    
    float vx2 = center.x - p1.x, vy2 = center.y - p1.y;
    float dot = (vx0 * vx1 + vy0 * vy1) * ((-vx0) * vx2 + (-vy0) * vy2);
    if(dot < 0 && !inside(p0) && !inside(p1)) {
      return null;
    }
    
//    println("vx", vx1, vy1, d0);
//    println(cross, h);
    
    float nx = -vy0 / d0, ny = vx0 / d0;
    int sign = cross < 0 ? 1 : -1;
    float hx = center.x + sign * h * nx;
    float hy = center.y + sign * h * ny;
    
//    println("n", nx, ny, sign, hx, hy);
    
    float d = sqrt(radius * radius - h * h);
    float tx = vx0 / d0, ty = vy0 / d0;
    float x = hx - tx * d;
    float y = hy - ty * d;
    
    if(pow(p1.x - x, 2) + pow(p1.y - y, 2) > d0 * d0) {
      x = hx + tx * d;
      y = hy + ty * d;
    }
    
//    println("a", d, tx, ty, x, y);
    
    return new Point(x, y);
  }
  public void move(float x, float y) {
    center.set(x, y);
  }  
}

float _lastUpdateTime;
boolean _isPaused = false;

ArrayList<Particle> _particles;
ArrayList<Point> _forces; // external force
ArrayList<Constraint> _constraints;
ArrayList<RigidBody> _bodies;
int _solverIterations = 4;
int _updateIterations = 2;

float _radius = 30;
Point _dragOrigin;

CircleBody _circle;
Point _crossed = new Point();

void setup() {
  size(500, 500);
  frameRate(30);
  
  _lastUpdateTime = millis();
  
  _particles = new ArrayList<Particle>();
  _forces = new ArrayList<Point>();
  _constraints = new ArrayList<Constraint>();
  _bodies = new ArrayList<RigidBody>();
  
//  Particle p = new Particle();
//  p.setPosition(width / 2, height / 2);
////  p.setVelocity(0, 13);
//  p.mass = 1;
//  _particles.add(p);
//  
//  Particle p2 = new Particle();
//  p2.setPosition(width * 3 / 4, height / 2);
////  p2.setVelocity(0, 13);
//  _particles.add(p2);
//  
//  float d = p.p0.dist(p2.p0);
//  DistanceConstraint dc = new DistanceConstraint(0, 1, d);
//  dc.k = 0.25f;
//  _constraints.add(dc);
//
//  Point p01 = new Point(width / 2, height / 2);
//  SpringConstraint sc = new SpringConstraint(0, p01);
//  sc.k = 1.0f;
//  _constraints.add(sc);
//  
//  Point p02 = new Point(width * 3 / 4, height * 3 / 4);
//  SpringConstraint sc2 = new SpringConstraint(1, p02);
//  sc2.k = 0.5f;
//  _constraints.add(sc2);

  _particles.add(new Particle(width / 2 - 150, height / 2));
  _particles.add(new Particle(width / 2 - 100, height / 2));
  _particles.add(new Particle(width / 2 - 50, height / 2));
  _particles.add(new Particle(width / 2, height / 2));
  
  _particles.add(new Particle(width / 2, height / 2 + 50));
  _particles.add(new Particle(width / 2 + 50, height / 2 + 50));
  _particles.add(new Particle(width / 2 + 100, height / 2 + 50));
  _particles.add(new Particle(width / 2 + 100, height / 2 + 100));
  _particles.add(new Particle(width / 2 + 100, height / 2 + 150));
  
  _particles.add(new Particle(width / 2 + 50, height / 2 + 150));
  _particles.add(new Particle(width / 2, height / 2 + 150));
  _particles.add(new Particle(width / 2 - 50, height / 2 + 150));
  _particles.add(new Particle(width / 2 - 100, height / 2 + 150));
  _particles.add(new Particle(width / 2 - 150, height / 2 + 150));
  
  _particles.add(new Particle(width / 2 - 150, height / 2 + 100));
  _particles.add(new Particle(width / 2 - 150, height / 2 + 50));
  
  for(int i = 0, len = _particles.size(); i < len; ++ i) {
    Particle pa = _particles.get(i);
    SpringConstraint sc = new SpringConstraint(i, pa.p0);
    _constraints.add(sc);
    //
    int j = (i + 1) % len;
    DistanceConstraint dc = new DistanceConstraint(i, j, 50);
    _constraints.add(dc);
  }
    
  //_forces.add(new Point(10, 0));
  
  _dragOrigin = new Point();
  
  _circle = new CircleBody(width / 2, height / 2, 30);
  _bodies.add(_circle);
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
  
  // finger
//  if(mousePressed) {
//    fill(255,255,255,64);
//  } else {
//    noFill();
//  }
//  stroke(255,255,255,128);
//  ellipse(mouseX, mouseY, 2 * _radius, 2 * _radius);

  fill(255,255,255,64);
  stroke(255);
  for(int i = 0, len = _bodies.size(); i < len; ++ i) {
    RigidBody b = _bodies.get(i);
    if(b instanceof CircleBody) {
      CircleBody c = (CircleBody) b;
      float r = 2 * c.radius;
      ellipse(c.center.x, c.center.y, r, r);
    }
  }
  
  fill(255,0,255);
  noStroke();
  ellipse(100, 100, 4, 4);
  ellipse(200, 100, 4, 4);
  noFill();
  stroke(255,0, 255);
  line(100, 100, 200, 100);
  if(_crossed != null) {
  fill(255,255,0);
  noStroke();
  ellipse(_crossed.x, _crossed.y, 8, 8);
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
//  Point[] ps = new Point[particles.size()];
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
  ArrayList<Constraint> constraints = generateCollisionConstraints(particles);//, ps);
  
  // loop solverIterations times
  //   projectConstraints(C1,...,CM+Mcoll ,p1,...,pN)
  // endloop
  for(int i = 0; i < _solverIterations; ++ i) {
    projectConstraints(_constraints, _particles);
    projectConstraints(constraints, _particles);
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

ArrayList<Constraint> generateCollisionConstraints(ArrayList<Particle> particles) {//, Point[] ps) {
  ArrayList<Constraint> res = new ArrayList<Constraint>();
  
  for(int i = 0, len = particles.size(); i < len; ++ i) {
    Particle pa = particles.get(i);
    Point p0 = pa.p0, p1 = pa.p1;
    for(int j = 0, len2 = _bodies.size(); j < len2; ++ j) {
      RigidBody b = _bodies.get(j);
      Point q = null;
      boolean s = b.inside(p0); // is static
      if(b.inside(p1)) {
        if(s) {
          q = b.closestPoint(p1);
        } else {
          q = b.intersection(p0, p1);
        } 
      } else {
        if(!s) {
          q = b.intersection(p0, p1);
        }
      }
      
      if(q == null) {
        continue;
      }
      
      Point n = b.normal(q);
      CollisionConstraint cc = new CollisionConstraint(i, q, n, s);
      res.add(cc);
    }
  }
  
  return res;
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
//  Particle pa = _particles.get(0);
//  pa.mass = 0;
//  pa.setPosition(mouseX, mouseY);

  _dragOrigin.set(mouseX, mouseY);
}

void mouseDragged() {
//  int i = _particles.size() - 1;
//  _particles.get(i).setPosition(mouseX, mouseY);
//  float dx = mouseX - width / 2, dy = mouseY - height / 2;
//  _forces.get(0).set(dx * 0.1, dy * 0.1);

//  Particle pa = _particles.get(0);
//  pa.setPosition(mouseX, mouseY);

  Point p0 = new Point(mouseX, mouseY);
//  PVector f = PVector.sub(p0, _dragOrigin);
//  f.normalize(); 
//  
//  for(int i = 0, len = _particles.size(); i < len; ++ i) {
//    Particle pa = _particles.get(i);
//    float d = pa.p0.dist(p0);
//    if(d > _radius) {
//      continue;
//    }
//    PVector n = PVector.div(PVector.sub(pa.p0, p0), d);
//    float v = _radius - d;
//    pa.setPosition(pa.p0.x + f.x * n.x * v, pa.p0.y + f.y * n.y * v);
//    pa.mass = 0; 
//  }

  _circle.move(p0.x, p0.y);
  
  Point p = _circle.intersection(new Point(100, 100), new Point(200, 100));
  println(p);
  if(p != null) {
    println(p.x, p.y);
  }
  
  _dragOrigin.set(p0);
}

void mouseReleased() {
//  Particle pa = _particles.get(0);
//  pa.mass = 1.0f;

  for(int i = 0, len = _particles.size(); i < len; ++ i) {
    if(_particles.get(i).mass == 0) {
      _particles.get(i).mass = 1;
    }
  }
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
