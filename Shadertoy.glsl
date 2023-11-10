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

struct Tore{
    float r;
    float R;
    int i;
};

struct Material{
    vec3 diffuse;
    vec3 specular;
    vec3 ambient;
    vec3 reflectionColor; // Couleur de réflexion
    float reflectionStrength; // Force de réflexion (0.0 pour aucune réflexion, 1.0 pour une réflexion totale)
};

Material normalMaterial(vec3 diff, vec3 spec, vec3 amb){
    Material mat;
    mat.diffuse = diff;
    mat.specular = spec;
    mat.ambient = amb;
    mat.reflectionColor = vec3(0.0);
    mat.reflectionStrength = 0.0;
    return mat;
}

Material reflectionMaterial(vec3 diff, vec3 spec, vec3 amb, vec3 reflColor, float reflStrength){
    Material mat;
    mat.diffuse = diff;
    mat.specular = spec;
    mat.ambient = amb;
    mat.reflectionColor = reflColor;
    mat.reflectionStrength = reflStrength;
    return mat;
}


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
    Tore tores[100];
    int nb_tores;
    vec3[100] lightPositions;
    int nb_lights;
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

void RotateCylinder(out Cylinder cyl, vec3 rot, float angle){
    // Translation du cylindre pour qu'il soit centré autour de l'origine
    vec3 center = (cyl.b + cyl.a) / 2.0;
    vec3 newa = cyl.a - center;
    vec3 newb = cyl.b - center;

    // Rotation dans l'espace local
    mat4 Rot = Rotation(rot, angle);
    vec3 newDirA = (Rot * vec4(newa, 0.0)).xyz;
    vec3 newDirB = (Rot * vec4(newb, 0.0)).xyz;

    // Retranslation pour ramener le cylindre à sa position d'origine
    newa = newDirA + center;
    newb = newDirB + center;

    cyl =  Cylinder(newa, newb, newDirA, cyl.r, cyl.i);
}

void RotateCapsule(out Capsule cap, vec3 rot, float angle){
    // Translation du cylindre pour qu'il soit centré autour de l'origine
    vec3 center = (cap.b + cap.a) / 2.0;
    vec3 newa = cap.a - center;
    vec3 newb = cap.b - center;

    // Rotation dans l'espace local
    mat4 Rot = Rotation(rot, angle);
    vec3 newDirA = (Rot * vec4(newa, 0.0)).xyz;
    vec3 newDirB = (Rot * vec4(newb, 0.0)).xyz;

    // Retranslation pour ramener le cylindre à sa position d'origine
    newa = newDirA + center;
    newb = newDirB + center;
    cap = Capsule(newa, newb, newDirB, cap.r, cap.i);
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

// Apply color model
// m : Material
// n : normal
vec3 Color(Material m,vec3 n){
    vec3 light=normalize(vec3(1,1,1));
 
    float diff=clamp(dot(n,light),0.,1.);
    vec3 col=m.diffuse*diff+vec3(.2,.2,.2);
    return col;
}

vec3 ColorUniform(Material m, vec3 n){
    return m.diffuse;
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

Material CheckerboardTexture(vec3 point) {
    // Définir la taille d'une case du damier
    float checkerSize = 1.0;

    // Ajouter une petite valeur d'epsilon aux coordonnées du point
    float epsilon = 0.1; // Ajustez cette valeur si nécessaire
    point += epsilon;

    // Calculer les indices de cases en x, y et z en utilisant des valeurs en virgule flottante
    float xIndex = floor(point.x / checkerSize);
    float yIndex = floor(point.y / checkerSize);
    float zIndex = floor(point.z / checkerSize);

    vec3 diff, spec, amb;
    amb = vec3(0.);

    // Alternance des cases diffuses et spéculaires en fonction des indices
    if ((int(xIndex) + int(yIndex) + int(zIndex))%2== 0) {
        // Case diffuse
        diff = vec3(0.8, 0.5, 0.3); 
        spec = vec3(0.1, 0.1, 0.1); 
    } else {
        // Case spéculaire
        diff = vec3(0.1, 0.1, 0.1);
        spec = vec3(0.2, 0.7, 0.8); 
    }

    return normalMaterial(diff, spec, amb);
}

// Compute color
// i : Texture index
// p : Point
Material Texture(vec3 p,int i){
    switch(i){
        case 0:{                                                        // CHECKBOARD
            return CheckerboardTexture(p);
        }
        case 1: return normalMaterial(vec3(.4,.0,.4),vec3(.2),vec3(.1));      // NORMAL VIOLET
        case 3: return normalMaterial(ColorConcentric(p),vec3(.0),vec3(.0));  // CONCENTRIC
        case 4: return normalMaterial(ColorRadial(p),vec3(.2),vec3(.1));      // RADIAL
        case 5: return normalMaterial(vec3(.8,.8,.8),vec3(.0),vec3(.0));      // UNIFORM
        case 6: return normalMaterial(vec3(.47,.15,.11),vec3(.1),vec3(.0));   // EARTH
        case 7: return normalMaterial(vec3(.2,.5,.3),vec3(.0),vec3(.0));      // UNIFORM GREEN
        case 8: return normalMaterial(vec3(.2,.2,.7),vec3(.0),vec3(.0));      // UNIFORM BLUE
        case 9: return normalMaterial(vec3(.8,.6,.2),vec3(.0),vec3(.0));      // UNIFORM ORANGE
        case 10: return reflectionMaterial(vec3(.8,.8,.8),vec3(.0),vec3(.0),vec3(.8,.8,.8),1.0); // REFLECTION
        default: return normalMaterial(vec3(0),vec3(0),vec3(0));
    }
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

bool IntersectEllipse(Ray ray, Ellipse e, vec3 rot, float angle, out Hit x) {
    mat4 Rot = Rotation(normalize(rot), angle);
    mat4 inv = inverse(Rot);
    vec3 rotCenter = (inv * vec4(e.c, 1.0)).xyz;
    vec3 rdd = (inv*vec4(ray.d,0.0)).xyz;
    vec3 roo = (inv*vec4(ray.o,1.0)).xyz;
    vec3 oc=roo-rotCenter;
    float a = dot(rdd/e.e, rdd/e.e);
    float b = 2.*dot(oc/e.e, rdd/e.e);
    float c = dot(oc/e.e, oc/e.e) - 1.;
    float d = b*b-4.*a*c;
    
    if(d>0.)
    {
        float t=min((-b-sqrt(d))/(2.*a), (-b+sqrt(d))/(2.*a));
        if(t>0.)
        {
            vec3 p=PointRot(roo,rdd,t);            
            x=Hit(t,normalize((Rot*vec4((p-rotCenter),0.0)).xyz), e.i);
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

bool IntersectBox(Ray ray, Box box, vec3 rot, float angle, out Hit x)
{
    mat4 Rot = Rotation(normalize(rot), angle);
    mat4 inv = inverse(Rot);
    vec3 center =(box.boxMax-box.boxMin)/2.+box.boxMin;
    vec3 rotCenter = (inv * vec4(center, 1.0)).xyz;
    
    vec3 Mm = (box.boxMax-box.boxMin)/2.;
    vec3 newMin = rotCenter-Mm;
    vec3 newMax = rotCenter+Mm;
    
    // ray to box
    vec3 rdd = (inv*vec4(ray.d,0.0)).xyz;
    vec3 roo = (inv*vec4(ray.o,1.0)).xyz;
    
    vec3 invDirection = 1.0 / rdd;
    vec3 tMin = (newMin - roo) * invDirection;
    vec3 tMax = (newMax - roo) * invDirection;
    vec3 t1 = min(tMin, tMax);
    vec3 t2 = max(tMin, tMax);
    float tNear = max(max(t1.x, t1.y), t1.z);
    float tFar = min(min(t2.x, t2.y), t2.z);
    if (tNear > tFar || tFar < 0.0) {
        return false; 
    }
    // Calcul normal
    vec4 normal = (tNear >0.0) ? vec4(tNear,(step(vec3(tNear),t1))) :
                                vec4(tFar,(step(t2,vec3(tFar))));
    normal.yzw = normalize((Rot*vec4(-sign(rdd)*normal.yzw,0.0)).xyz);
    x = Hit(tNear, normal.yzw, 1);
    return true;
}

bool IntersectTore(Ray ray, Tore tor, vec3 rot,vec3 trs, float angle, out Hit x){
    mat4 rt = Rotation(normalize(rot), angle);
    mat4 tr = Translate(trs);
    mat4 txi = tr * rt;
    mat4 txx = inverse(txi);
    vec3 rdd = (txx*vec4(ray.d,0.0)).xyz;
    vec3 roo = (txx*vec4(ray.o,1.0)).xyz;

    // vec3 oc = ray.o - tor.ct;
    // float po = 1.0;
    float Ra2 = tor.r*tor.r;
    float ra2 = tor.R*tor.R;
    // float m = dot(ray.o,ray.o);
    float m = dot(roo,roo);
    // float n = dot(ray.o,ray.d);
    float n = dot(roo,rdd);

    float k = (m - ra2 - Ra2)/2.0;
    float k3 = n;
  
    
    float k2 = n*n + Ra2*rdd.z*rdd.z + k;
    float k1 = k*n + Ra2*roo.z*rdd.z;
    float k0 = k*k + Ra2*roo.z*roo.z - Ra2*ra2;
    
    float c2 = 2.0*k2 - 3.0*k3*k3;
    float c1 = k3*(k3*k3 - k2) + k1;
    float c0 = k3*(k3*(-3.0*k3*k3 + 4.0*k2) - 8.0*k1) + 4.0*k0;
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;
    float Q = c2*c2 + c0;
    float R = 3.0*c0*c2 - c2*c2*c2 - c1*c1;
    float h = R*R - Q*Q*Q;
    float z = 0.0;
    if(h<0.0){
        float sQ = sqrt(Q);
        z = 2.0*sQ*cos(acos(R/(sQ*Q)) / 3.0);
    } else {
        float sQ = pow( sqrt(h) + abs(R), 1.0/3.0 );
        z = sign(R)*abs( sQ + Q/sQ );
    }
    z = c2 - z;
    float d1 = z - 3.0*c2;
    float d2 = z*z - 3.0*c0;
    if( abs(d1) < 1.0e-4){
        if(d2<0.0){
            return false;
        }
        d2 = sqrt(d2);
    } else {
        if(d1 < 0.0){
            return false;
        }
        d1 = sqrt( d1/2.0 );
        d2 = c1/d1;
    }
    float r = 1e20;
    h = d1*d1 - z + d2;
    if(h>0.0){
        h = sqrt(h);
        float t1 = -d1 - h - k3;
        // t1 = (po<0.0)?2.0/t1:t1;
        float t2 = -d1 + h - k3;
        // t2 = (po<0.0)?2.0/t2:t2;
        if(t1>0.0){
            r=t1;
        }
        if(t2>0.0){
            r=min(r,t2);
        }
    }
    h = d1*d1 - z - d2;
    if(h>0.0){
        h=sqrt(h);
        float t1 = d1 - h - k3;
        float t2 = d1 + h - k3;
        if(t1>0.0){
            r=min(r,t1);
        }
        if(t2>0.0){
            r=min(r,t2);
        }
    }
    // vec3 p = Point(ray,r);
    vec3 p = roo + r*rdd;
    vec3 nor = p*(dot(p,p) - tor.R*tor.R - tor.r*tor.r*vec3(1.0,1.0,-1.0));
    vec4 tmp = vec4(r, nor);
    tmp.yzw = (txi*vec4(tmp.yzw,0.0)).xyz;
    x = Hit(tmp.x, normalize(tmp.yzw), tor.i);
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

void addToreToScene(out Scene sc, Tore tor){
    sc.tores[sc.nb_tores] = tor;
    sc.nb_tores++;
}

void addLightToScene(out Scene sc, vec3 position){
    sc.lightPositions[sc.nb_lights] = position;
    sc.nb_lights++;
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
        if(IntersectEllipse(ray,scene.ellipses[i],vec3(0.0,0.0,1.0),iTime, current)&&current.t<x.t){
            x = current;
            ret = true;
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
        if(IntersectBox(ray,scene.boxes[i],vec3(0.0,0.0,1.0),iTime, current)&&current.t<x.t){
            x = current;
            ret = true;
        }
}

void IntersectTore(Ray ray, out Hit x, out bool ret, out Hit current){
    for(int i = 0; i < scene.nb_tores; i++)
        if(IntersectTore(ray,scene.tores[i],vec3(-1.0,0.0,0.0), vec3(-2.,-2.,0.), iTime,current)&&current.t<x.t){
        x=current;
        ret=true;
    }
}

void RotateScene(out Scene sc, vec3 rot, float angle){
    for(int i = 0; i < sc.nb_cylinders; i++)
        RotateCylinder(sc.cylinders[i], rot, angle);

    for(int i = 0; i < sc.nb_capsules; i++)
        RotateCapsule(sc.capsules[i], rot, angle);
}

// Scene intersection
// ray : The ray
//   x : Returned intersection information
bool Intersect(Ray ray,out Hit x){
    // ROTATION
    const vec3 ROT = vec3(0.0,0.0,1.0);
    float angle = iTime;
    
    x=Hit(1000.,vec3(0),-1);
    Hit current;
    bool ret=false;

    IntersectSpheres(ray, x, ret, current);
    IntersectPlanes(ray, x, ret, current);
    IntersectEllipses(ray, x, ret, current);
    IntersectTruncatedCylinders(ray, x, ret, current);
    IntersectCapsules(ray, x, ret, current);
    IntersectBoxes(ray, x, ret, current);
    IntersectTore(ray, x, ret, current);
    
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

bool IsInShadow(vec3 p, vec3 lightPosition, Scene scene)
{
    vec3 lightDir = normalize(lightPosition - p);
    Ray shadowRay = Ray(p, lightDir);
    Hit x;
    bool idx = Intersect(shadowRay, x);
    return idx;
}

vec3 CalculatePhongLighting(vec3 point, Material material, vec3 normal, Scene sc) {
    vec3 ambient = material.ambient; 
    vec3 diffuse = material.diffuse;
    vec3 specular = material.specular; 

    vec3 totalLighting = vec3(0.0);
    float epsilon = 0.5; // Ajouter un petit décalage

    for (int i = 0; i < sc.nb_lights; i++) {
        // Calculate lighting components (diffuse and specular) for each light source
        vec3 lightDir = normalize(sc.lightPositions[i] - point);

        // Ajouter epsilon au point pour corriger la précision numérique
        vec3 pointWithEpsilon = point + epsilon * lightDir;

        if(IsInShadow(pointWithEpsilon, sc.lightPositions[i], sc)) {
            continue; // Ignorer la source lumineuse si elle est en ombre
        }

        // Calcul de la distance entre le point et la source lumineuse
        float distance = length(sc.lightPositions[i] - point);

        // Calcul de l'atténuation en fonction de la distance (loi de l'inverse du carré)
        float coeffLight = 0.25;
        float attenuation = 1.0 / pow(distance,coeffLight);

        float diff = max(dot(normal, lightDir), 0.0);

        // Appliquer l'atténuation à la composante diffuse
        diff *= attenuation;

        // Calculate the reflect direction for the specular component
        vec3 reflectDir = reflect(-lightDir, normal);

        // Appliquer l'atténuation à la composante spéculaire
        float specularAttenuation = attenuation; // Vous pouvez ajuster cela si nécessaire

        // Combine ambient, diffuse, and specular for the current light source
        vec3 lighting = ambient + diffuse * diff + specular * specularAttenuation;

        // Accumulate the lighting from all light sources
        totalLighting += lighting;
    }
    return totalLighting;
}

vec3 CalculateReflection(vec3 incidentDir, vec3 normal, vec3 intersectionPoint, Hit intersectionInfo, Scene scene, int maxDepth) {
    if (maxDepth <= 0) {
        return vec3(0.0); // Arrêtez les réflexions après un certain nombre d'itérations
    }

    // Calculer la direction réfléchie
    vec3 reflectedDir = reflect(incidentDir, normal);

    // Créer un rayon réfléchi à partir du point d'intersection
    Ray reflectedRay;
    reflectedRay.o = intersectionPoint;
    reflectedRay.d = reflectedDir;

    // Recherche d'intersection avec les objets environnants
    Hit reflectionHit;
    bool hasReflectionIntersection = Intersect(reflectedRay, reflectionHit);

    // S'il y a une intersection, calculez la couleur de la réflexion récursivement
    if (hasReflectionIntersection) {
        Material reflectionMaterial = Texture(intersectionPoint, reflectionHit.i);
        // Calculez l'éclairage pour le matériau réfléchissant (vous pouvez utiliser CalculatePhongLighting ou autre)
        vec3 reflectionLighting = CalculatePhongLighting(intersectionPoint, reflectionMaterial, reflectionHit.n, scene);

        // Récursivement calculez la couleur de la réflexion
        vec3 reflectionColor = reflectionLighting /** CalculateReflection(reflectedDir, reflectionHit.n, intersectionPoint, reflectionHit, scene, maxDepth - 1)*/;

        return reflectionColor;
    } else {
        return vec3(0.0); // Pas d'intersection, couleur de fond
    }
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
    
        // Calculate Phong lighting
        vec3 lighting = CalculatePhongLighting(p, mat, x.n, scene);

        // For reflectionMaterial
        if (mat.reflectionStrength > 0.0) {
            vec3 reflectionColor = CalculateReflection(ray.d, x.n, p, x, scene, 1);
            lighting = mix(lighting, reflectionColor, mat.reflectionStrength);
        }

        return lighting;
    }
    else
        return Background(ray.d);
    return vec3(0);
}

void initializeScene(out Scene sc){
    addSphereToScene(sc, Sphere(vec3(5.,-3.,1.5),1.,6));
    addSphereToScene(sc, Sphere(vec3(0.,0.,4.),1.,10));
    addPlaneToScene(sc, Plane(vec3(0.,0.,1.),vec3(0.,0.,-2.),0));
    addEllipseToScene(sc, Ellipse(vec3(-5.,2.,3.0), vec3(2., 3., 1.), 8));
    addCylinderToScene(sc, Cylinder(vec3(4.,-2.,2), vec3(2.,-2.,2), vec3(1.,0.,0.), vec3(0.,1.,0.),3));
    addCylinderToScene(sc, Cylinder(vec3(4.,1.,4), vec3(4.,1.,2), vec3(0.,0.,1.), vec3(0.,1.,0.),5));
    addCapsuleToScene(sc, Capsule(vec3(-1.,3.,1), vec3(-1.,5.,1), vec3(0.,1.,0.), 1.,7));
    addBoxToScene(sc, Box(vec3(-4.,-4.,4.),vec3(-2.,-2.,2.),9));
    addToreToScene(sc, Tore(1., 0.5, 9));

    //addLightToScene(sc, vec3(10.,0.,10.));
    addLightToScene(sc, vec3(0.,5.,7.));
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
    vec3 color=Shade(Ray(ro,rd));
    fragColor = vec4(color, 1.0);
}
