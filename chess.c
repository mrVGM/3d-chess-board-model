#include <math.h>
#include <GL/glfw.h>
#include <stdio.h>

typedef enum{false, true} bool;

bool running() {
    return( !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam( GLFW_OPENED) );
}

void init() {
    int width, height;

    glfwInit();
    if( !glfwOpenWindow( 640, 480, 0, 0, 0, 0, 0, 0, GLFW_WINDOW ) ) return;

    glfwGetWindowSize( &width, &height );
    height = height > 0 ? height : 1;

    glViewport( 0, 0, width, height );
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 65.0f, (GLfloat)width/(GLfloat)height, 1.0f, 100.0f );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt(0.0f, -10.0f, 6.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f );
}
void finit(){
    glfwTerminate();
}
typedef struct Point {
    float x;
    float y;
    float z;
} Point;

float vLen(Point p1) {
    return sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
}

Point vectProd(Point p1, Point p2) {
    return (Point) {p1.y*p2.z - p1.z*p2.y, -p1.x*p2.z + p1.z*p2.x, p1.x*p2.y - p1.y*p2.x};
}

Point vect(Point p1, Point p2) {
    return (Point) {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
}
Point midPoint(Point p1, Point p2, float t) {
    Point v = vect(p1, p2);
    return (Point) {p1.x + v.x*t, p1.y + v.y*t, p1.z + v.z*t};
}

Point norm(Point p) {
    float len = vLen(p);
    return (Point) {p.x/len, p.y/len, p.z/len};
}

void makeTrMat(Point a, float* res) {
    Point xAx = norm(a);
    Point yAx = (Point) {0,1,0};
    Point zAx = vectProd(xAx, yAx);
    res[0] = xAx.x;
    res[1] = yAx.x;
    res[2] = zAx.x;

    res[3] = xAx.y;
    res[4] = yAx.y;
    res[5] = zAx.y;

    res[6] = xAx.z;
    res[7] = yAx.z;
    res[8] = zAx.z;
}

Point trByMat(float* m, Point p) {
    Point res;
    res.x = m[0]*p.x + m[1]*p.y + m[2]*p.z;
    res.y = m[3]*p.x + m[4]*p.y + m[5]*p.z;
    res.z = m[6]*p.x + m[7]*p.y + m[8]*p.z;
    return res;
}

void covr(Point a, Point b, Point c, Point d) {
    Point base1[30];
    Point base2[30];

    Point nv1[30];
    Point nv2[30];

    float mat[9];
    makeTrMat(vect(a, b), mat);
    float r = vLen(vect(a ,b))/2;
    Point center = midPoint(a, b, 0.5);
    int i;
    for (i = 0; i < 30; i++) {
        Point tmp = (Point) {cos(2*M_PI*i/30.0)*r, sin(2*M_PI*i/30.0)*r, 0};
        tmp = trByMat(mat, tmp);
        base1[i] = (Point) {center.x + tmp.x, center.y + tmp.y, center.z + tmp.z};
        nv1[i] = norm(tmp);
    }
    makeTrMat(vect(c, d), mat);
    r = vLen(vect(c, d))/2;
    center = midPoint(c, d, 0.5);
    for (i = 0; i < 30; i++) {
        Point tmp = (Point) {cos(2*M_PI*i/30.0)*r, sin(2*M_PI*i/30.0)*r, 0};
        tmp = trByMat(mat, tmp);
        base2[i] = (Point) {center.x + tmp.x, center.y + tmp.y, center.z + tmp.z};
        nv2[i] = norm(tmp);
    }
    Point c1 = vect(a, b);
    Point c2 = vect(c, d);

    glBegin(GL_TRIANGLES);
    for (i = 0; i < 29; i++) {
        glNormal3f(nv1[i].x, nv1[i].y, nv1[i].z);
        glVertex3f(base1[i].x, base1[i].y, base1[i].z);

        glNormal3f(nv1[i+1].x, nv1[i+1].y, nv1[i+1].z);
        glVertex3f(base1[i+1].x, base1[i+1].y, base1[i+1].z);

        glNormal3f(nv2[i+1].x, nv2[i+1].y, nv2[i+1].z);
        glVertex3f(base2[i+1].x, base2[i+1].y, base2[i+1].z);


        glNormal3f(nv1[i].x, nv1[i].y, nv1[i].z);
        glVertex3f(base1[i].x, base1[i].y, base1[i].z);

        glNormal3f(nv2[i].x, nv2[i].y, nv2[i].z);
        glVertex3f(base2[i].x, base2[i].y, base2[i].z);

        glNormal3f(nv2[i+1].x, nv2[i+1].y, nv2[i+1].z);
        glVertex3f(base2[i+1].x, base2[i+1].y, base2[i+1].z);
    }

    glNormal3f(nv1[29].x, nv1[29].y, nv1[29].z);
    glVertex3f(base1[29].x, base1[29].y, base1[29].z);

    glNormal3f(nv1[0].x, nv1[0].y, nv1[0].z);
    glVertex3f(base1[0].x, base1[0].y, base1[0].z);

    glNormal3f(nv2[0].x, nv2[0].y, nv2[0].z);
    glVertex3f(base2[0].x, base2[0].y, base2[0].z);

    glNormal3f(nv1[29].x, nv1[29].y, nv1[29].z);
    glVertex3f(base1[29].x, base1[29].y, base1[29].z);

    glNormal3f(nv2[29].x, nv2[29].y, nv2[29].z);
    glVertex3f(base2[29].x, base2[29].y, base2[29].z);

    glNormal3f(nv2[0].x, nv2[0].y, nv2[0].z);
    glVertex3f(base2[0].x, base2[0].y, base2[0].z);

    glEnd();
}

void bezHelper(Point* p, int *n, float t) {
    int i;
    int h = *n;
    for (i = 0; i < h-1; i++) {
        Point op;
        op = (Point) {(1-t)*p[i].x + t*p[i+1].x, (1-t)*p[i].y + t*p[i+1].y, (1-t)*p[i].z + t*p[i+1].z};
        p[i] = op;
    }
    *n = *n-1;
}
Point bez(Point* p, int n, float t) {
    Point* tmp = malloc(n*sizeof(Point));
    memcpy(tmp, p, n*sizeof(Point));
    while (n != 1) {
        bezHelper(tmp, &n, t);
    }
    Point res = tmp[0];
    free (tmp);
    return res;
}
void drBez(Point* p, int n, int pts, Point* dest) {
    int i;
    for (i = 0; i < pts; i++) {
        Point tmp = bez(p, n, (float)i/(pts-1));
        dest[i] = (Point) {tmp.x, tmp.y, tmp.z};
    }
}

void drHalfSphere(float ang, float r, float h, Point pos) {
    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);

    glPushMatrix();
    glRotatef(ang, 0, 1, 0);

    glPushMatrix();
    glScalef(r,r,h);

    glBegin(GL_TRIANGLES);

    int i;
    for (i = 0; i < 15; i++) {
        int j;
        for (j = 0; j < 30; j++) {
            glNormal3f(cos(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));
            glVertex3f(cos(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));

            glNormal3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));
            glVertex3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));

            glNormal3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));
            glVertex3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));

            glNormal3f(cos(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));
            glVertex3f(cos(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*i/15.0), sin(M_PI*i/15.0));

            glNormal3f(cos(2*M_PI*j/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));
            glVertex3f(cos(2*M_PI*j/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*j/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));

            glNormal3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));
            glVertex3f(cos(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(2*M_PI*(j+1)/30.0)*cos(M_PI*(i+1)/15.0), sin(M_PI*(i+1)/15.0));
        }
    }
    glEnd();

    glPopMatrix();
    glPopMatrix();
    glPopMatrix();
}


void drawUnitCylinder ( )
{
    int na = 30;
    float dalpha = 2*M_PI/na;
    float dx = cos(dalpha);
    float dy = sin(dalpha);

    float x1 = 1;
    float y1 = 0;

    int i;
    for(i=0; i<na; i++)
    {

        float x2 = x1*dx-y1*dy;
        float y2 = y1*dx+x1*dy;

        // Ðèñóâàíå íà îêîëíà ñòåíà
        glBegin( GL_POLYGON );
            glNormal3f( x1, y1, 0 );
            glVertex3f( x1, y1, 0 );
            glVertex3f( x1, y1, 1 );
            glNormal3f( x2, y2, 0 );
            glVertex3f( x2, y2, 1 );
            glVertex3f( x2, y2, 0 );
        glEnd();

        // Ðèñóâàíå íà ïàð÷å îò äîëíàòà îñíîâà
        glBegin( GL_POLYGON );
            glNormal3f(  0,  0, -1 );
            glVertex3f(  0,  0,  0 );
            glVertex3f( x1, y1,  0 );
            glVertex3f( x2, y2,  0 );
        glEnd();

        // Ðèñóâàíå íà ïàð÷å îò ãîðíàòà îñíîâà
        glBegin( GL_POLYGON );
            glNormal3f(  0,  0, 1 );
            glVertex3f(  0,  0, 1 );
            glVertex3f( x1, y1, 1 );
            glVertex3f( x2, y2, 1 );
        glEnd();

        x1 = x2;
        y1 = y2;
    }
}

void drCone (float r)
{
    int na = 30;
    float dalpha = 2*M_PI/na;
    float dx = cos(dalpha);
    float dy = sin(dalpha);

    float x1 = 1;
    float y1 = 0;

    int i;
    for(i=0; i<na; i++)
    {

        float x2 = x1*dx-y1*dy;
        float y2 = y1*dx+x1*dy;

        // Ðèñóâàíå íà îêîëíà ñòåíà
        glBegin( GL_POLYGON );
            glNormal3f( x1, y1, 0 );
            glVertex3f( x1, y1, 0 );
            glVertex3f( r*x1, r*y1, 1 );
            glNormal3f( x2, y2, 0 );
            glVertex3f( r*x2, r*y2, 1 );
            glVertex3f( x2, y2, 0 );
        glEnd();

        // Ðèñóâàíå íà ïàð÷å îò äîëíàòà îñíîâà
        glBegin( GL_POLYGON );
            glNormal3f(  0,  0, -1 );
            glVertex3f(  0,  0,  0 );
            glVertex3f( x1, y1,  0 );
            glVertex3f( x2, y2,  0 );
        glEnd();

        // Ðèñóâàíå íà ïàð÷å îò ãîðíàòà îñíîâà
        glBegin( GL_POLYGON );
            glNormal3f(  0,  0, 1 );
            glVertex3f(  0,  0, 1 );
            glVertex3f( r*x1, r*y1, 1 );
            glVertex3f( r*x2, r*y2, 1 );
        glEnd();

        x1 = x2;
        y1 = y2;
    }
}


Point spherical( float alpha, float beta, float r )
{
    Point p;
    p.x = r*cos(alpha)*cos(beta);
    p.y = r*sin(alpha)*cos(beta);
    p.z = r*sin(beta);
    return p;
}

void drawSmoothSphere()
{
    float x = 0;
    float y = 0;
    float z = 0;
    float r = 1;
    int na = 15;
    int nb = 15;
    float alpha;
    float dalpha = 2.0*M_PI/na;
    float beta;
    float dbeta = 1.0*M_PI/nb;

    beta = M_PI/2;
    int j;
    for(j=0; j<nb; j++, beta-=dbeta)
    {
        alpha = 0;
        int i;
        for(i=0; i<na; i++, alpha+=dalpha)
        {
            Point p;
            glBegin( GL_POLYGON );

            p = spherical(alpha,beta,1);
            glNormal3f(p.x,p.y,p.z);
            glVertex3f(x+r*p.x,y+r*p.y,z+r*p.z);

            p = spherical(alpha+dalpha,beta,1);
            glNormal3f(p.x,p.y,p.z);
            glVertex3f(x+r*p.x,y+r*p.y,z+r*p.z);

            p = spherical(alpha+dalpha,beta-dbeta,1);
            glNormal3f(p.x,p.y,p.z);
            glVertex3f(x+r*p.x,y+r*p.y,z+r*p.z);

            p = spherical(alpha,beta-dbeta,1);
            glNormal3f(p.x,p.y,p.z);
            glVertex3f(x+r*p.x,y+r*p.y,z+r*p.z);

            glEnd( );
        }
    }
}

void drTetr() {
    glBegin(GL_TRIANGLES);

    glNormal3f(0,0,-1);
    glVertex3f(-1, -1, 0);
    glVertex3f(1, -1, 0);
    glVertex3f(1, 1, 0);
    glVertex3f(-1, -1, 0);
    glVertex3f(-1, 1, 0);
    glVertex3f(1, 1, 0);

    Point nv1 = norm(vectProd((Point) {1, 0, 0}, (Point){1,1,1}));
    Point nv2 = norm(vectProd((Point) {0,1,0}, (Point) {-1, 1, 1}));

    glNormal3f(nv1.x, nv1.y, nv1.z);
    glVertex3f(-1,-1,0);
    glVertex3f(1,-1,0);
    glVertex3f(0,0,1);

    glNormal3f(nv2.x, nv2.y, nv2.z);
    glVertex3f(1,-1,0);
    glVertex3f(1,1,0);
    glVertex3f(0,0,1);


    glNormal3f(-nv1.x, -nv1.y, nv1.z);
    glVertex3f(1,1,0);
    glVertex3f(-1,1,0);
    glVertex3f(0,0,1);

    glNormal3f(-nv2.x, -nv2.y, nv2.z);
    glVertex3f(-1,-1,0);
    glVertex3f(-1,1,0);
    glVertex3f(0,0,1);

    glEnd();
}

void mane(Point *b, Point *p, int n, float offset) {
    int i;
    int koef = 1;
    if (offset < 0) {
        koef = -1;
    };



    glBegin(GL_TRIANGLES);
    for (i = 0; i < n-1; i++) {
        Point p1 = (Point) {b[i].x, b[i].y + offset, b[i].z};
        Point p2 = (Point) {b[i+1].x, b[i+1].y + offset, b[i+1].z};
        Point p3 = (Point) {p[i+1].x, p[i+1].y, p[i+1].z};


        Point nv = norm(vectProd(vect(p3, p2), vect(p3, p1)));


        glNormal3f(nv.x, koef*nv.y, nv.z);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);


        p1 = (Point) {b[i].x, b[i].y + offset, b[i].z};
        p2 = (Point) {p[i].x, p[i].y, p[i].z};
        p3 = (Point) {p[i+1].x, p[i+1].y, p[i+1].z};

        nv = norm(vectProd(vect(p1, p2), vect(p1, p3)));

        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);

    }
    glEnd();
}

void drMane() {
    Point p[13], b[10];
    p[	0	]=(Point){	0.06	,	0	,	-0.16	};
    p[	1	]=(Point){	0.16	,	0	,	-0.22	};
    p[	2	]=(Point){	0.29	,	0	,	-0.51	};
    p[	3	]=(Point){	0.23	,	0	,	-0.72	};
    p[	4	]=(Point){	0.17	,	0	,	-0.93	};
    p[	5	]=(Point){	0.3	,	0	,	-0.91	};
    p[	6	]=(Point){	0.26	,	0	,	-1.05	};
    p[	7	]=(Point){	0.16	,	0	,	-1.41	};
    p[	8	]=(Point){	0	,	0	,	-1.48	};
    p[	9	]=(Point){	-0.14	,	0	,	-1.79	};
    p[	10	]=(Point){	-0.27	,	0	,	-1.68	};
    p[	11	]=(Point){	-0.52	,	0	,	-1.66	};
    p[	12	]=(Point){	-0.65	,	0	,	-1.51	};

    b[	0	]=(Point){	-0.11	,	0	,	-0.19	};
    b[	1	]=(Point){	-0.02	,	0	,	-0.29	};
    b[	2	]=(Point){	0.04	,	0	,	-0.47	};
    b[	3	]=(Point){	0.03	,	0	,	-0.6	};
    b[	4	]=(Point){	0.03	,	0	,	-0.73	};
    b[	5	]=(Point){	0.02	,	0	,	-1	};
    b[	6	]=(Point){	-0.09	,	0	,	-1.1	};
    b[	7	]=(Point){	-0.21	,	0	,	-1.2	};
    b[	8	]=(Point){	-0.2	,	0	,	-1.41	};
    b[	9	]=(Point){	-0.61	,	0	,	-1.4	};

    Point x[30], y[30];

    drBez(p, 13, 30, x);
    drBez(b, 10, 30, y);
    int i;


    mane(y, x, 29, 0.1);
    mane(y, x, 29, -0.1);

    glBegin(GL_TRIANGLES);
    Point p1 = (Point) {y[28].x, y[28].y + 0.1, y[28].z};
    Point p2 = (Point) {y[28].x, y[28].y - 0.1, y[28].z};
    Point p3 = (Point) {x[28].x, x[28].y, x[28].z};
    Point nor = norm(vectProd(vect(p3, p1), vect(p3, p2)));
    glNormal3f(-nor.x, -nor.y, -nor.z);
    glVertex3f(p1.x, p1.y, p1.z);
    glVertex3f(p2.x, p2.y, p2.z);
    glVertex3f(p3.x, p3.y, p3.z);
    glEnd();
}

void knight() {
    glPushMatrix();
    glTranslatef(0,0,0.2);
    glRotatef(180, 1,0,0);
    glPushMatrix();
    Point p[7], v[4];
    p[	0	]=(Point){	-0.01	,	0	,	0.37	};
    p[	1	]=(Point){	-0.22	,	0	,	0.1	};
    p[	2	]=(Point){	0.61	,	0	,	-1.06	};
    p[	3	]=(Point){	-0.42	,	0	,	-1.47	};
    p[	4	]=(Point){	-0.68	,	0	,	-1.57	};
    p[	5	]=(Point){	-1.02	,	0	,	-1.22	};
    p[	6	]=(Point){	-1.06	,	0	,	-1.17	};

    v[	0	]=(Point){	-0.97	,	0	,	0.38	};
    v[	1	]=(Point){	-1.38	,	0	,	-0.25	};
    v[	2	]=(Point){	0.32	,	0	,	-1.14	};
    v[	3	]=(Point){	-0.84	,	0	,	-0.87	};


    Point curve1[30];
    drBez(p, 7, 30, curve1);
    Point curve2[30];
    drBez(v, 4, 30, curve2);

    Point centerPoint = vect(midPoint(curve1[0], curve2[0], 0.5), (Point) {0,0,0});
    glTranslatef(centerPoint.x, centerPoint.y, centerPoint.z);

    int i;
    for (i = 0; i < 29; i++) {
        covr(curve2[i], curve1[i], curve2[i+1], curve1[i+1]);
    }
    Point tmp = norm(vect(curve1[29], curve2[29]));
    float angle = acos(tmp.x);
    if (tmp.z > 0) {
        angle = -1*angle;
    }
    angle = 180*angle/M_PI;
    drHalfSphere(angle, vLen(vect(curve1[29], curve2[29]))/2, 0.3, midPoint(curve1[29], curve2[29], 0.5));

    glPushMatrix();
    glTranslatef(-0.1,0,0.2);
    drMane();
    glPopMatrix();
    glPopMatrix();
    glPopMatrix();
    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-0.5,-0.15,1.7);
    glScalef(0.05, 0.05, 0.05) ;
    drawSmoothSphere(1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-0.5,0.15,1.7);
    glScalef(0.05, 0.05, 0.05) ;
    drawSmoothSphere(1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-0.25,-0.1,1.78);
    glRotatef(20, 1,0,0);
    glScalef(0.05,0.05,0.25);

    drTetr();

    glPopMatrix();

    glPushMatrix();
    glTranslatef(-0.25,0.1,1.78);
    glRotatef(-20, 1,0,0);
    glScalef(0.05,0.05,0.25);

    drTetr();

    glPopMatrix();
}



void drawOxyz()
{
    // Ðèñóâàíå íà ëú÷èòå
    glBegin( GL_LINES );
        // OX
        glVertex3f( 0.0, 0.0, 0.0 );
        glVertex3f( 5.0, 0.0, 0.0 );
        // OY
        glVertex3f( 0.0, 5.0, 0.0 );
        glVertex3f( 0.0, 0.0, 0.0 );
        // OZ
        glVertex3f( 0.0, 0.0, 5.0 );
        glVertex3f( 0.0, 0.0, 0.0 );
    glEnd();

    // Ðèñóâàíå íà ñòðåëêèòå
    glBegin( GL_TRIANGLES );
        // OX
        glVertex3f( 5.0, 0.0, 0.0 );
        glVertex3f( 4.5, 0.2, 0.0 );
        glVertex3f( 4.5,-0.2, 0.0 );
        // OY
        glVertex3f( 0.0, 5.0, 0.0 );
        glVertex3f( 0.2, 4.5, 0.0 );
        glVertex3f(-0.2, 4.5, 0.0 );
        // OZ
        glVertex3f( 0.0,  0.0,  5.0 );
        glVertex3f(+0.14,-0.14, 4.5 );
        glVertex3f(-0.14,+0.14, 4.5 );
    glEnd();
}

void pawn() {
    glPushMatrix();
    glScalef(1,1,0.9);

    glPushMatrix();
    glTranslatef(0,0,0.2);
    glPushMatrix();
    glScalef(0.5,0.5,0.9);
    drCone(0.3);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0,0,1.2);
    glScalef(0.35,0.35,0.35);
    drawSmoothSphere();
    glPopMatrix();
    glPopMatrix();

    glPopMatrix();

    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();


}

void bishop() {
    glPushMatrix();
    glScalef(1,1,0.9);

    glPushMatrix();
    glTranslatef(0,0,0.2);

    glPushMatrix();
    glScalef(1,1,1.25);


    glPushMatrix();
    glScalef(0.5,0.5,0.9);
    drCone(0.3);
    glPopMatrix();


    glPushMatrix();
    glTranslatef(0,0,1.25);
    glScalef(0.35,0.35,0.35);
    drawSmoothSphere();
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0,0,0.87);
    glScalef(0.35,0.35,0.05);
    drawUnitCylinder();
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.55);
    glScalef(0.1,0.1,0.7);
    drCone(0.01);
    glPopMatrix();
    glPopMatrix();


    glPopMatrix();


    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();

}

void drawSolidCube (float x,float y,float z,float a)
{
    a = a/2;
    glBegin( GL_POLYGON );
        glNormal3f(-1.0, 0.0, 0.0);
        glVertex3f(x-a, y-a, z-a);
        glVertex3f(x-a, y-a, z+a);
        glVertex3f(x-a, y+a, z+a);
        glVertex3f(x-a, y+a, z-a);
    glEnd();
    glBegin( GL_POLYGON );
        glNormal3f(+1.0, 0.0, 0.0);
        glVertex3f(x+a, y-a, z-a);
        glVertex3f(x+a, y-a, z+a);
        glVertex3f(x+a, y+a, z+a);
        glVertex3f(x+a, y+a, z-a);
    glEnd();
    glBegin( GL_POLYGON );
        glNormal3f(0.0, -1.0, 0.0);
        glVertex3f(x-a, y-a, z-a);
        glVertex3f(x-a, y-a, z+a);
        glVertex3f(x+a, y-a, z+a);
        glVertex3f(x+a, y-a, z-a);
    glEnd();
    glBegin( GL_POLYGON );
        glNormal3f(0.0, +1.0, 0.0);
        glVertex3f(x-a, y+a, z-a);
        glVertex3f(x-a, y+a, z+a);
        glVertex3f(x+a, y+a, z+a);
        glVertex3f(x+a, y+a, z-a);
    glEnd();
    glBegin( GL_POLYGON );
        glNormal3f(0.0, 0.0, -1.0);
        glVertex3f(x-a, y-a, z-a);
        glVertex3f(x-a, y+a, z-a);
        glVertex3f(x+a, y+a, z-a);
        glVertex3f(x+a, y-a, z-a);
    glEnd();
    glBegin( GL_POLYGON );
        glNormal3f(0.0, 0.0, +1.0);
        glVertex3f(x-a, y-a, z+a);
        glVertex3f(x-a, y+a, z+a);
        glVertex3f(x+a, y+a, z+a);
        glVertex3f(x+a, y-a, z+a);
    glEnd();
}

void unitCube() {
    drawSolidCube(0,0,0.5,1);
}

void fortress() {
    glPushMatrix();
    glScalef(5,5,1);
    unitCube();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-2,-2,1);
    glPushMatrix();
    unitCube();
    glTranslatef(2,0,0);
    unitCube();
    glTranslatef(2,0,0);
    unitCube();

    glTranslatef(0,2,0);
    unitCube();
    glTranslatef(-4,0,0);
    unitCube();

    glTranslatef(0,2,0);
    unitCube();
    glTranslatef(2,0,0);
    unitCube();
    glTranslatef(2,0,0);
    unitCube();
    glPopMatrix();
    glPopMatrix();
}

void rook() {
    glPushMatrix();
    glTranslatef(0,0,0.2);
    glPushMatrix();
    glScalef(0.5,0.5,1.4);
    drCone(0.53);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(0,0,1.4);
    glScalef(0.13,0.13,0.13);
    fortress();
    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();

}

void queen() {
    glPushMatrix();
    glScalef(1,1,1.1);

    glPushMatrix();
    glTranslatef(0,0,0.2);


    glPushMatrix();
    glScalef(0.5,0.5,1.5);
    drCone(0.53);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.55);
    glScalef(0.265, 0.265, 0.3);
    drCone(1.5);

    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.5);
    glScalef(0.35,0.35,0.04);
    drawUnitCylinder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.93);
    glScalef(0.1, 0.1, 0.1);
    drawSmoothSphere();
    glPopMatrix();

    glPopMatrix();

    glPopMatrix();

    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();
}

void cross() {
    glPushMatrix();
    glScalef(1,1,4);
    drawSolidCube(0, 0, 0.5, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,2);
    glScalef(3,1,1);
    drawSolidCube(0, 0, 0.5, 1);
    glPopMatrix();
}

void king() {
    glPushMatrix();
    glScalef(1,1,1.1);

    glPushMatrix();
    glTranslatef(0,0,0.2);


    glPushMatrix();
    glScalef(0.5,0.5,1.5);
    drCone(0.53);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.55);
    glScalef(0.265, 0.265, 0.3);
    drCone(1.5);

    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.5);
    glScalef(0.35,0.35,0.04);
    drawUnitCylinder();
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,1.75);
    glScalef(0.1, 0.1, 0.1);
    cross();
    glPopMatrix();

    glPopMatrix();
    glPopMatrix();

    glPushMatrix();
    glScalef(0.6,0.6,0.2);
    drawUnitCylinder();
    glPopMatrix();
}

int colorN;
void chCol() {
    if (colorN%2 == 0) {
        glColor3ub(50,50,50);
    }
    else {
        glColor3ub(200,200,200);
    }
    colorN++;
}

void board() {
    glPushMatrix();
    glTranslatef(-4.9, -4.9, 0);
    glPushMatrix();
    glTranslatef(0,0,-1);
    glScalef(1.4, 1.4, 1);
    glPushMatrix();
    glTranslatef(3.5,3.5,-0.1);
    glScalef(8.5,8.5,1);
    glColor3ub(255, 255, 0);
    unitCube();
    glPopMatrix();
    int i;
    colorN = 0;
    chCol();
    glPushMatrix();
    for (i = 0; i < 8; i++) {
        int j;
        for (j = 0; j < 8; j++) {
            unitCube();
            chCol();
            glTranslatef(1,0,0);
        }
        glTranslatef(-8, 1, 0);
        chCol();
    }
    glPopMatrix();
    glPopMatrix();

    glColor3ub(255,255,255);
    glPushMatrix();
    rook();
    glTranslatef(1.4, 0, 0);

    glPushMatrix();
    glRotatef(-90,0,0,1);
    knight();
    glPopMatrix();

    glTranslatef(1.4, 0, 0);
    bishop();

    glTranslatef(1.4, 0, 0);
    queen();

    glTranslatef(1.4, 0, 0);
    king();

    glTranslatef(1.4, 0, 0);
    bishop();

    glTranslatef(1.4, 0, 0);
    glPushMatrix();
    glRotatef(-90,0,0,1);
    knight();
    glPopMatrix();

    glTranslatef(1.4, 0, 0);
    rook();

    glTranslatef(-9.8, 1.4, 0);
    for (i = 0; i < 8; i++) {
        pawn();
        glTranslatef(1.4, 0,0);
    }

    glTranslatef(-11.2, 7, 0);

    glColor3ub(20,20,20);
    for (i = 0; i < 8; i++) {
        pawn();
        glTranslatef(1.4, 0,0);
    }

    glTranslatef(-11.2, 1.4, 0);


    rook();
    glTranslatef(1.4, 0, 0);

    glPushMatrix();
    glRotatef(90,0,0,1);
    knight();
    glPopMatrix();

    glTranslatef(1.4, 0, 0);
    bishop();

    glTranslatef(1.4, 0, 0);
    queen();

    glTranslatef(1.4, 0, 0);
    king();

    glTranslatef(1.4, 0, 0);
    bishop();

    glTranslatef(1.4, 0, 0);
    glPushMatrix();
    glRotatef(90,0,0,1);
    knight();
    glPopMatrix();

    glTranslatef(1.4, 0, 0);
    rook();
    glPopMatrix();
    glPopMatrix();
}

int main()
{
    init();

    glEnable( GL_DEPTH_TEST );
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);

    glEnable(GL_SMOOTH);
	glShadeModel(GL_SMOOTH);

    glEnable(GL_NORMALIZE);

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );

    while( running() )
    {
        glClear( GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT );
        glRotatef( 0.05, 0.4, -0.2, 0.7);

        board();

        glfwSwapBuffers();
    }

    return 0;
}

