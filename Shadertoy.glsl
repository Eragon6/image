struct Sphere{
    vec3 c;// Center
    float r;// Radius
    int i;// Texture Id
};

struct Plane{
    vec3 n;// Normal
    vec3 p;// Point
    int i;// Texture Id
};

struct Hit{
    float t;// Intersection depth
    vec3 n;// Normal
    int i;// Texture Id
};

struct Ray{
    vec3 o;// Origin
    vec3 d;// Direction
};

struct Tore{
    vec3 Ro;
    float r;
    float R;
    int i;
};

struct Box{
    vec3 boxMin;
    vec3 boxMax;
    int i;
};

struct Capsule{
    vec3 a;
    vec3 b;
    vec3 dir;
    float r;
    int i;
};


struct Cylinder{
    vec3 a;
    vec3 b;
    vec3 dir;
    vec3 r;
    int i;
};

struct Ellipse{
    vec3 c;
    vec3 e; // e(a demi grand axe, b demi petit, c demi distance focale)
    int i;
};

struct Material{
    vec3 d;// Diffuse
};


struct Scene{
    Sphere spheres[100];
    int nb_spheres;
    Plane planes[100];
    int nb_planes;
    Ellipse ellipses[100];
    int nb_ellipses;
    Cylinder cylinders[100];
    int nb_cylinders;
    Capsule capsules[100];
    int nb_capsules;
    Box boxes[100];
    int nb_boxes;
};

Scene scene;

mat4 Translate(vec3 tr){
    mat4 translate_matrix = mat4( 1.0   , 0.0   , 0.0   , 0.0,
                                  0.0   , 1.0   , 0.0   , 0.0,
                                  0.0   , 0.0   , 1.0   , 0.0,
                                  tr.x   , tr.y   , tr.z   , 1.0 );
    return translate_matrix;
}

mat4 Rotation(vec3 rot, float theta){
    float s = sin(theta);
    float c = cos(theta);
    float compc = 1.-cos(theta);
    
    mat4 rotation_matrix = mat4(rot.x*rot.x*compc+c,       rot.x*rot.y*compc-s*rot.z, rot.x*rot.z*compc+s*rot.y, 0.0,
                                rot.x*rot.y*compc+s*rot.z, rot.y*rot.y*compc+c,       rot.y*rot.z*compc-s*rot.x, 0.0,
                                rot.x*rot.z*compc-s*rot.y, rot.y*rot.z*compc+s*rot.x, rot.z*rot.z*compc+c,       0.0,
                                0.0,                       0.0,                       0.0,                       1.0);
    return rotation_matrix;
}

Cylinder RotateCylinder(Cylinder cyl, vec3 rot, float angle){
    // Translation du cylindre pour qu'il soit centré autour de l'origine
    vec3 center = (cyl.b - cyl.a) / 2.0 + cyl.a;
    vec3 newa = cyl.a - center;
    vec3 newb = cyl.b - center;

    // Rotation dans l'espace local
    mat4 Rot = Rotation(rot, angle);
    vec3 newDirA = (Rot * vec4(newa, 0.0)).xyz;
    vec3 newDirB = (Rot * vec4(newb, 0.0)).xyz;

    // Retranslation pour ramener le cylindre à sa position d'origine
    newa = newDirA + center;
    newb = newDirB + center;

    return Cylinder(newa, newb, newDirA, cyl.r, cyl.i);
}

Capsule RotateCapsule(Capsule cap, vec3 rot, float angle){
    // Translation du cylindre pour qu'il soit centré autour de l'origine
    vec3 center = (cap.b - cap.a) / 2.0 + cap.a;
    vec3 newa = cap.a - center;
    vec3 newb = cap.b - center;

    // Rotation dans l'espace local
    mat4 Rot = Rotation(rot, angle);
    vec3 newDirA = (Rot * vec4(newa, 0.0)).xyz;
    vec3 newDirB = (Rot * vec4(newb, 0.0)).xyz;

    // Retranslation pour ramener le cylindre à sa position d'origine
    newa = newDirA + center;
    newb = newDirB + center;
    return Capsule(newa, newb, newDirB, cap.r, cap.i);
}

// TODO : Rotation d'ellipse
Ellipse RotateEllipse(Ellipse ell, vec3 rot, float angle){
    return ell;
}

Box RotateBox(Box box, vec3 rot, float angle){
    return box;
}

float Checkers(in vec2 p){
    // Filter kernel
    vec2 w=fwidth(p)+.001;
    // Box box filter
    vec2 i=2.*(abs(fract((p-.5*w)*.5)-.5)-abs(fract((p+.5*w)*.5)-.5))/w;
    // xor pattern
    return.5-.5*i.x*i.y;
}

// Hemisphere direction
// seed : Integer seed, from 0 to N
//    n : Direction of the hemisphere
vec3 Hemisphere(int seed,vec3 n){
    float a=fract(sin(176.19*float(seed)));// Uniform randoms
    float b=fract(sin(164.19*float(seed)));
    
    float u=2.*3.1415*a;// Random angle
    float v=acos(2.*b-1.);// Arccosine distribution to compensate at poles
    
    vec3 d=vec3(cos(u)*cos(v),sin(u)*cos(v),sin(v));// Direction
    if(dot(d,n)<0.){d=-d;}// Hemisphere
    
    return d;
}

// Compute point on ray
// ray : The ray
//   t : Distance
vec3 Point(Ray ray,float t){
    return ray.o+t*ray.d;
}

vec3 PointRot(vec3 ro, vec3 rd,float t){
    return ro+t*rd;
}

// Compute color
// i : Texture index
// p : Point
Material Texture(vec3 p,int i){
    if(i==1)
    {
        return Material(vec3(.8,.5,.4));
    }
    else if(i==0)
    {
        // compute checkboard
        float f=Checkers(.5*p.xy);
        vec3 col=vec3(.4,.5,.7)+f*vec3(.1);
        return Material(col);
    }
    else if(i==5)
    {
        return Material(vec3(.8,.8,.4));
    }
    return Material(vec3(0));
}

// Sphere intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectSphere(Ray ray,Sphere sph, out Hit x){
    vec3 oc=ray.o-sph.c;
    float b=dot(oc,ray.d);
    float c=dot(oc,oc)-sph.r*sph.r;
    float d=b*b-c;
    if(d>0.)
    {
        float t=-b-sqrt(d);
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            x=Hit(t,normalize(p-sph.c),sph.i);
            
            return true;
        }
    }
    return false;
}

// Plane intersection
// ray : The ray
//   x : Returned intersection information
bool IntersectPlane(Ray ray,Plane pl,out Hit x){
    float t=-dot(ray.o-pl.p,pl.n)/dot(ray.d,pl.n);
    if(t>0.)
    {
        
        x=Hit(t,vec3(0,0,1),pl.i);
        return true;
    }
    return false;
}


bool IntersectEllipse(Ray ray, Ellipse e,  out Hit x) {
    vec3 oc=ray.o-e.c;
    float a = dot(ray.d/e.e, ray.d/e.e);
    float b = 2.*dot(oc/e.e, ray.d/e.e);
    float c = dot(oc/e.e, oc/e.e) - 1.;
    float d = b*b-4.*a*c;
    
    if(d>0.)
    {
        float t=min((-b-sqrt(d))/(2.*a), (-b+sqrt(d))/(2.*a));
        if(t>0.)
        {
            vec3 p=Point(ray,t);
            x=Hit(t,normalize((p-e.c)/dot(e.e, e.e)), e.i);
            return true;
        }
    }
    return false;    
}

bool IntersectTruncatedCylinder(Ray ray, Cylinder cyl, out Hit x){
    float d2 = dot(ray.d, ray.d);
    vec3 u = normalize(cyl.dir);
    vec3 oa = ray.o - cyl.a;
    float a = d2 - dot(dot(ray.d, u), dot(ray.d, u));
    float b = 2.0 * (dot(oa, ray.d) - dot(oa, u) * dot(ray.d, u));
    float c = dot(oa, oa) - (dot(dot(oa, u), dot(oa, u))) - (dot(cyl.r, cyl.r));

    float delta = b * b - 4.0 * a * c;
    if (delta > 0.0) {
        float t1 = (-b - sqrt(delta)) / (2.0 * a);
        float t2 = (-b + sqrt(delta)) / (2.0 * a);

        float t = min(t1, t2);

        // Vérifier si l'intersection est à l'intérieur du cylindre tronqué
        vec3 p = Point(ray, t);
        float height = dot(p - cyl.b, u);
        if (t > 0.0 && height >= 0.0 && height <= length(cyl.b-cyl.a)) {
            x = Hit(t, normalize(p - cyl.a - dot(p - cyl.a, u) * u), cyl.i);
            return true;
        }
    }

    // Intersections avec les disques tronqués aux extrémités du cylindre tronqué
    float tTop = (dot(cyl.b - ray.o, u)) / dot(ray.d, u);
    float tBottom = (dot(cyl.a - ray.o, u)) / dot(ray.d, u);

    vec3 pTop = Point(ray, tTop);
    vec3 pBottom = Point(ray, tBottom);

    bool hitTop = false;
    bool hitBottom = false;

    if (tTop > 0.0 && length(pTop - cyl.b) <= length(cyl.r)) {
        x = Hit(tTop, normalize(cyl.b - cyl.a), cyl.i);
        hitTop = true;
        if (tBottom > 0.0 && length(pBottom - cyl.a) <= length(cyl.r) && tTop > tBottom) {
            x = Hit(tBottom, normalize(cyl.a - cyl.b), cyl.i);
        }
        return true;
    }
    if (tBottom > 0.0 && length(pBottom - cyl.a) <= length(cyl.r)) {
        x = Hit(tBottom, normalize(cyl.a - cyl.b), cyl.i);
        return true;
    }
    return false;
}

bool IntersectCapsule(Ray ray, Capsule cap, out Hit x){
    Hit hitTop, hitBottom;
    bool hitTopSphere, hitBottomSphere;
    
    float d2 = dot(ray.d, ray.d);
    vec3 u = normalize(cap.dir);
    vec3 oa = ray.o - cap.a;
    float a = d2 - dot(dot(ray.d, u), dot(ray.d, u));
    float b = 2.0 * (dot(oa, ray.d) - dot(oa, u) * dot(ray.d, u));
    float c = dot(oa, oa) - (dot(dot(oa, u), dot(oa, u))) - (cap.r * cap.r);
    float delta = b * b - 4.0 * a * c;
    if (delta > 0.0) {
        float t1 = (-b - sqrt(delta)) / (2.0 * a);
        float t2 = (-b + sqrt(delta)) / (2.0 * a);

        float t = min(t1, t2);

        // Vérifier si l'intersection est à l'intérieur du cylindre tronqué
        vec3 p = Point(ray, t);
        float height = dot(p - cap.a, u);
        if (t > 0.0 && height >= 0.0 && height <= length(cap.b-cap.a)) {
            x = Hit(t, normalize(p - cap.a - dot(p - cap.a, u) * u), cap.i);
            return true;
        }
        vec3 captb = cap.a;
        bool comp = dot(p,cap.dir) < dot(captb,cap.dir);
        if(height>length(cap.b-cap.a)){
            oa = ray.o-cap.b;
            captb = cap.b;
            comp = dot(p,cap.dir) > dot(captb,cap.dir);
        }
        
        b=dot(oa, ray.d);
        c=dot(oa,oa)-cap.r*cap.r;
        delta=b*b-c;
        if (delta >0.){
            float t = -b-sqrt(delta);
            if (t > 0.)
            {
                vec3 p = Point(ray, t);
                if (comp){
                    x=Hit(t,normalize(p-captb),cap.i);
                    return true;
                }
            }  
        }
    }
    return false;
}

bool IntersectBox(Ray ray, Box box, out Hit x)
{
    vec3 invDirection = 1.0 / ray.d;
    vec3 tMin = (box.boxMin - ray.o) * invDirection;
    vec3 tMax = (box.boxMax - ray.o) * invDirection;

    vec3 t1 = min(tMin, tMax);
    vec3 t2 = max(tMin, tMax);

    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);

    if (tNear > tFar || tFar < 0.0) {
        return false; // Pas d'intersection
    }

    vec3 intersectionPoint = Point(ray, tNear);
    vec3 normal;

    // Déterminer la normale de la face de la boîte touchée
    if (tNear == t1.x) normal = normalize(intersectionPoint-(box.boxMin+box.boxMax)/2.);
    else if (tNear == t1.y) normal = normalize(intersectionPoint-(box.boxMin+box.boxMax)/2.);
    else if (tNear == t1.z) normal = normalize(intersectionPoint-(box.boxMin+box.boxMax)/2.);
    /*else if (tNear == t2.x) normal = normalize(intersectionPoint-(box.boxMin+box.boxMax)/2.);
    else if (tNear == t2.y) normal = normalize(intersectionPoint-(box.boxMin-box.boxMax)/2.);
    else normal = normalize(intersectionPoint-(box.boxMin+box.boxMax)/2.);*/

    x = Hit(tNear, normal, 1);
    return true;
}

void addSphereToScene(out Scene sc, Sphere sph){
    sc.spheres[sc.nb_spheres] = sph;
    sc.nb_spheres++;
}

void addPlaneToScene(out Scene sc, Plane pl){
    sc.planes[sc.nb_planes] = pl;
    sc.nb_planes++;
}

void addEllipseToScene(out Scene sc, Ellipse ell){
    sc.ellipses[sc.nb_ellipses] = ell;
    sc.nb_ellipses++;
}

void addCylinderToScene(out Scene sc, Cylinder cyl){
    sc.cylinders[sc.nb_cylinders] = cyl;
    sc.nb_cylinders++;
}

void addCapsuleToScene(out Scene sc, Capsule cap){
    sc.capsules[sc.nb_capsules] = cap;
    sc.nb_capsules++;
}

void addBoxToScene(out Scene sc, Box box){
    sc.boxes[sc.nb_boxes] = box;
    sc.nb_boxes++;
}

void IntersectSpheres(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_spheres; i++)
        if(IntersectSphere(ray,scene.spheres[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void IntersectPlanes(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_planes; i++)
        if(IntersectPlane(ray,scene.planes[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void IntersectEllipses(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_ellipses; i++)
        if(IntersectEllipse(ray,scene.ellipses[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void IntersectTruncatedCylinders(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_cylinders; i++)
        if(IntersectTruncatedCylinder(ray,scene.cylinders[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void IntersectCapsules(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_capsules; i++)
        if(IntersectCapsule(ray,scene.capsules[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void IntersectBoxes(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_boxes; i++)
        if(IntersectBox(ray,scene.boxes[i],current)&&current.t<x.t){
            x=current;
            ret=true;
        }
}

void RotateScene(out Scene sc, vec3 rot, float angle){
    for(int i = 0; i < sc.nb_cylinders; i++)
        sc.cylinders[i] = RotateCylinder(sc.cylinders[i], rot, angle);

    for(int i = 0; i < sc.nb_capsules; i++)
        sc.capsules[i] = RotateCapsule(sc.capsules[i], rot, angle);

    for(int i = 0; i < sc.nb_ellipses; i++)
        sc.ellipses[i] = RotateEllipse(sc.ellipses[i], rot, angle);

    for(int i = 0; i < sc.nb_boxes; i++)
        sc.boxes[i] = RotateBox(sc.boxes[i], rot, angle);
}

// Scene intersection
// ray : The ray
//   x : Returned intersection information
bool Intersect(Ray ray,out Hit x){
    // ROTATION
    vec3 ROT = vec3(0.0,0.0,1.0);
    float ANGLE = iTime;
    // Cylinder cyl = RotationCylinder(scene.cylinders[0],ROT,ANGLE);
    // Capsule cap = RotationCapsule(scene.capsules[0],ROT,ANGLE);
    // ell = RotationEllipse(ell,ROT,ANGLE);
    
    x=Hit(1000.,vec3(0),-1);
    Hit current;
    bool ret=false;

    IntersectSpheres(ray, x, ret, current);
    IntersectPlanes(ray, x, ret, current);
    IntersectEllipses(ray, x, ret, current);
    IntersectTruncatedCylinders(ray, x, ret, current);
    IntersectCapsules(ray, x, ret, current);
    IntersectBoxes(ray, x, ret, current);
    
    return ret;
}

vec3 Background(vec3 rd){
    return mix(vec3(.8,.8,.9),vec3(.7,.7,.8),rd.z);
}

// Camera rotation matrix
// ro : Camera origin
// ta : Target point
mat3 setCamera(in vec3 ro,in vec3 ta){
    vec3 cw=normalize(ta-ro);
    vec3 cp=vec3(0,0,1);
    vec3 cu=-normalize(cross(cw,cp));
    vec3 cv=-normalize(cross(cu,cw));
    return mat3(cu,cv,cw);
}

// Apply color model
// m : Material
// n : normal
vec3 Color(Material m,vec3 n){
    vec3 light=normalize(vec3(1,1,1));
 
    float diff=clamp(dot(n,light),0.,1.);
    vec3 col=m.d*diff+vec3(.2,.2,.2);
    return col;
}

vec3 ColorUniform(Material m, vec3 n){
    return m.d;
}


vec3 ColorCheckers(Material mat, vec3 p, vec3 n){
    // Size of a checker cell
    float CHECKER_SIZE = 1.0;

    // calculating the value using Checkers
    float checkerValue = Checkers(p.xy / CHECKER_SIZE);

    vec3 color1 = vec3(0.9, 0.9, 0.9);
    vec3 color2 = vec3(0.3, 0.3, 0.3);

    // Choose the color using checkerValue
    vec3 checkerColor = mix(color1, color2, checkerValue);

    // Add a lighting using the normal
    vec3 light=normalize(vec3(1,1,1));
    float diff=clamp(dot(n,light),0.,1.);
    // apply the lighting
    vec3 finalColor = checkerColor * (diff + 0.2);

    return finalColor;
}

vec3 HSVtoRGB(vec3 c){
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 ColorConcentric(vec3 p){
   // Centre de la texture
    const vec3 CENTER = vec3(0.0, 0.0, 0.);

    // Distance from the center
    float distance = length(p - CENTER);

    // Variation with the distance
    float hue = mod(distance, 1.0); 
    vec3 color = HSVtoRGB(vec3(hue, 1.0, 1.0)); 
    return color;
}

vec3 ColorRadial(vec3 p){
    const vec3 AXIS = normalize(vec3(0.0, 0.0, 1.0));

    // Distance from the axis
    float distanceToAxis = length(cross(p, AXIS));

    // Color variation
    float hue = distanceToAxis / 2.0;

    float sat = 1.0;
    float lum = 1.0; 
    
    vec3 rgbColor = HSVtoRGB(vec3(hue, sat, lum));
    
    return rgbColor;
}

// Rendering
vec3 Shade(Ray ray){
    // Intersect contains all the geo detection
    Hit x;
    bool idx=Intersect(ray,x);
    
    if(idx)
    {
        vec3 p=Point(ray,x.t);
        Material mat=Texture(p,x.i);
        
        switch(x.i)
        { 
            case 2: { return ColorCheckers(mat,p,x.n); }
            case 3: { return ColorConcentric(p); }
            case 4: { return ColorRadial(p); }
            case 5: { return ColorUniform(mat,x.n); }
            default: { return Color(mat,x.n); }
        }
    }
    else
    {
        return Background(ray.d);
    }
    
    return vec3(0);
}

void initializeScene(out Scene sc){
    addSphereToScene(sc, Sphere(vec3(0.,0.,1.),1.,3));
    addSphereToScene(sc, Sphere(vec3(1.,2.,0.),1.,4));
    addPlaneToScene(sc, Plane(vec3(0.,0.,1.),vec3(0.,0.,0.),2));
    addEllipseToScene(sc, Ellipse(vec3(-5.,2.,3.0), vec3(2., 3., 1.), 2));
    addCylinderToScene(sc, Cylinder(vec3(4.,-2.,2), vec3(2.,-2.,2), vec3(1.,0.,0.), vec3(0.,1.,0.),3));
    addCylinderToScene(sc, Cylinder(vec3(4.,1.,4), vec3(4.,1.,2), vec3(0.,0.,1.), vec3(0.,1.,0.),5));
    addCapsuleToScene(sc, Capsule(vec3(-1.,3.,1), vec3(-1.,5.,1), vec3(0.,1.,0.), 1.,4));
    addBoxToScene(sc, Box(vec3(-4.,-4.,4.),vec3(-2.,-2.,2.), 1));
}


void mainImage(out vec4 fragColor,in vec2 fragCoord){
    // From uv which are the pixel coordinates in [0,1], change to [-1,1] and apply aspect ratio
    vec2 uv=(-iResolution.xy+2.*fragCoord.xy)/iResolution.y;
    
    // Mouse control
    vec2 mouse=iMouse.xy/iResolution.xy;
    
    // Ray origin
    vec3 ro=12.*normalize(vec3(sin(2.*3.14*mouse.x),cos(2.*3.14*mouse.x),1.4*(mouse.y-.1)));
    vec3 ta=vec3(0.,0.,1.5);
    mat3 ca=setCamera(ro,ta);
    
    // Ray
    vec3 rd=ca*normalize(vec3(uv.xy*tan(radians(22.5)),1.));

    initializeScene(scene);
    RotateScene(scene, vec3(0.0,0.0,1.0), iTime);
    
    // Render
    vec3 col=Shade(Ray(ro,rd));
    
    fragColor=vec4(col,1.);
}
