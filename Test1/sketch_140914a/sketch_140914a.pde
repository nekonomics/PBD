class Point extends PVector {
  public Point() {
    this(0, 0);
  }
  public Point(float x, float y) {
    super(x, y);
  }
}

class Particle {
  public Point p;
  public Point v;
  public float mass = 1.0;
  public Particle() {
    p = new Point();
    v = new Point();
  }
  public void set(float x, float y, float vx, float vy) {
    setPosition(x, y);
    setVelocity(vx, vy);
  }
  public void setPosition(float x, float y) {
    p.set(x, y);
  }
  public void setVelocity(float vx, float vy) {
    v.set(vx, vy);
  }
  public float x() {
    return p.x;
  }
  public float y() {
    return p.y;
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
  public Point[] positions;
  
  // stiffness
  // 0...1
  public float k = 1.0f;
  
  // equality or inequality
  //   equality: eval() == 0
  // inequality: eval() >= 0
  public boolean isEquality = true;
  
  // evaluation function
  // R^{3n} -> R
  public float eval() {
    return 0;
  }
  
  public void setPosition(int i, Point p) {
    positions[i] = p;
  }
  
  public Point getPosition(int i) {
    return positions[i];
  }
}

class DistanceConstraint extends Constraint {
  public float distance; 
  public DistanceConstraint(int i1, int i2, float distance) {
    this.distance = distance;
    n = 2;
    indices = new int[] { i1, i2 };
    positions = new Point[n];
  }  
  public float eval() {
    Point p1 = positions[0];
    Point p2 = positions[1];
//    Point p = new Point(p1.x - p2.x, p1.y - p2.y);
    float d = p1.dist(p2);
    float c = d - distance;
    
    if(!isEquality && c >= 0) {
      return c;
    }
    
    float w1 = 1.0 / 1.0;
    float w2 = 1.0 / 1.0;
    float s = c / (w1 + w2);
    PVector n = PVector.div(PVector.sub(p1, p2), d);
    
    Point dp1 = new Point(-s * w1 * n.x, -s * w1 * n.y);
    Point dp2 = new Point( s * w2 * n.x,  s * w2 * n.y);
    
    float dt = 1.0 / 30; // TODO
    positions[0].set(p1.x + dt * dp1.x, p1.y + dt * dp1.y);
    positions[1].set(p2.x + dt * dp2.x, p2.y + dt * dp2.y);
    
    return c;
  }
}

float _lastUpdateTime;

ArrayList<Particle> _particles;
ArrayList<Point> _forces; // external force
ArrayList<Constraint> _constraints;
int _solveIterations = 2;

void setup() {
  size(500, 500);
  frameRate(30);
  
  _lastUpdateTime = millis();
  
  _particles = new ArrayList<Particle>();
  _forces = new ArrayList<Point>();
  _constraints = new ArrayList<Constraint>();
  
  Particle p = new Particle();
  p.setPosition(width / 2, height / 2);
  p.setVelocity(0, 13);
  _particles.add(p);
  
  Particle p2 = new Particle();
  p2.setPosition(p.x() + 100, p.y());
  p2.setVelocity(0, 13);
  _particles.add(p2);
  
  DistanceConstraint dc = new DistanceConstraint(0, 1, 50);
  _constraints.add(dc);
  
  _forces.add(new Point(10, 0));
}

void draw() {
  background(0);
  
  float time = millis();
  float dt = (time - _lastUpdateTime) / 1000.0f;
  _lastUpdateTime = time;
  
  Point force = new Point(0, 0);
  for(int i = 0; i < _forces.size(); ++ i) {
    Point f = _forces.get(i);
    force.x += f.x;
    force.y += f.y;
  }
  
//  for(int i = 0; i < _particles.size(); ++ i) {
//    Particle p = _particles.get(i);
//    updateParticle(p, force, dt);
//  }

  updateParticles(_particles, force, dt);
  
  for(int i = 0; i < _particles.size(); ++ i) {
    Particle p = _particles.get(i);
    drawParticle(p);
    if(i == 0) {
      text("p=(" + p.p + ") \nv=(" + p.v + ")", 10, 10);
    }
  }
}

void updateParticles(ArrayList<Particle> particles, Point force, float dt) {
  // forall vertices i do vi ← vi + ∆t wi fext (xi)
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    float vx = p.vx() + dt * force.x / p.mass;
    float vy = p.vy() + dt * force.y / p.mass;
    p.setVelocity(vx, vy);
  }
  // dampVelocities(v1 , . . . , vN )

  // forall vertices i do pi ← xi + ∆t vi
  Point[] ps = new Point[particles.size()];
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    float x = p.x() + dt * p.vx();
    float y = p.y() + dt * p.vy();
    ps[i] = new Point(x, y);
  }
  // forall vertices i do generateCollisionConstraints(xi → pi)
  
  // loop solverIterations times
  //   projectConstraints(C1,...,CM+Mcoll ,p1,...,pN)
  // endloop
  for(int i = 0; i < _solveIterations; ++ i) {
    projectConstraints(_constraints, ps);
  }
  
  // forall vertices i
  //   vi ←(pi−xi)/∆t
  //   xi ← pi
  // endfor
  for(int i = 0; i < particles.size(); ++ i) {
    Particle p = particles.get(i);
    float x = ps[i].x;
    float y = ps[i].y;
    float vx = (x - p.x()) / dt;
    float vy = (y - p.y()) / dt;
    p.set(x, y, vx, vy);
  }
  
  // velocityUpdate(v1 , . . . , vN )
  
}

ArrayList<Constraint> generateCollisionConstraints() {
  return null;
}

void projectConstraints(ArrayList<Constraint> constraints, Point[] positions) {
  for(int i = 0; i < constraints.size(); ++ i) {
    Constraint c = constraints.get(i);
    for(int j = 0; j < c.n; ++ j) {
      c.setPosition(j, positions[c.indices[j]]);
    }
    float r = c.eval();
    for(int j = 0; j < c.n; ++ j) { 
      positions[c.indices[j]].set(c.getPosition(j));      
    }
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

void drawParticle(Particle p) {
  float x = p.x(), y = p.y();
  fill(255);
  noStroke();
  ellipse(x, y, 8, 8);
  noFill();
  stroke(255,128,0);
  line(x, y, x + p.vx(), y + p.vy());
}
