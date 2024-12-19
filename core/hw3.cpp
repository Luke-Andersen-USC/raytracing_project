/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Luke Andersen
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <cmath> 
#include <algorithm>
#include <iostream>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = 1;

#define WIDTH 640
#define HEIGHT 480
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vector3 {
    float x, y, z;

    Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}
    Vector3(const double input[3]) : x(input[0]), y(input[1]), z(input[2]) {}
    Vector3(const Vector3& other) : x(other.x), y(other.y), z(other.z) {}

    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    Vector3 operator-() const {
        return Vector3(-x, -y, -z);
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(const float other) const {
        return Vector3(x * other, y * other, z * other);
    }

   Vector3 cross(const Vector3& other) const {
      return Vector3(
          y * other.z - z * other.y,
          z * other.x - x * other.z,
          x * other.y - y * other.x
      );
    }

    float dot(const Vector3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    void normalize() {
        float mag = sqrt(x * x + y * y + z * z);
        if (mag > 1e-8) {
            x /= mag;
            y /= mag;
            z /= mag;
        }
        else
        {
          x = y = z = 0.0f;
        }
    }

    float dist(const Vector3& other) const {
        return sqrt(pow(other.x - x, 2.0) + pow(other.y - y, 2.0) + pow(other.z - z, 2.0));
    }

    float length() const {
      return sqrt(x * x + y * y + z * z);
    }
};


class Object {
public:
  virtual ~Object() {}
  virtual Vector3 GetCenter() const { return Vector3(0.0f, 0.0f, 0.0f); };

};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

class Triangle : public Object
{

public:
  Vertex v[3];
  Vector3 GetCenter() const
  {
    float x = (v[0].position[0] +v[1].position[0] + v[2].position[0]) / 3.0f;
    float y = (v[0].position[1] +v[1].position[1] + v[2].position[1]) / 3.0f;
    float z = (v[0].position[2] +v[1].position[2] + v[2].position[2]) / 3.0f;
    return Vector3(x, y, z);
  }
};

class Sphere : public Object
{
public:
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;

  Vector3 GetCenter() const
  {
    return Vector3(position);
  }
};

struct Light
{
  double position[3];
  double color[3];
};

struct Ray
{
  Vector3 origin;
  Vector3 direction;

  Ray(Vector3 inOrigin = Vector3(0.0f, 0.0f, 0.0f), Vector3 dest = Vector3(0.0f, 0.0f, 0.0f))
  {
    origin = inOrigin;
    direction = dest - origin;
    direction.normalize();
  }
};

struct Intersection
{
  Object* object;

  Vector3 position;
  Vector3 toCamera;
  Vector3 normal;
  Vector3 color_diffuse;
  Vector3 color_specular;
  double shininess;
};


Ray rays[HEIGHT][WIDTH];
unsigned char redDisplay[HEIGHT][WIDTH];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

void draw_scene()
{
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      unsigned char r = buffer[y][x][0]; // modify
      unsigned char g = buffer[y][x][1]; // modify
      unsigned char b = buffer[y][x][2]; // modify
      plot_pixel(x, y, r, g, b);
    }
    glEnd();
    glFlush();
  }
  printf("Ray tracing completed.\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else 
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

Intersection CalculateSphereIntersection(const Ray& ray, Sphere& sphere);
Intersection  CalculateTriangleIntersection(const Ray& ray, Triangle& triangle);

bool CheckInShadow(const Intersection& intersection, Light light, bool shouldPrintStuff = false)
{
  Vector3 lightPosition(light.position);
  Ray shadowRay(intersection.position, lightPosition);
  shadowRay.origin = shadowRay.origin + shadowRay.direction * 0.0005f;

  for(int i = 0; i < num_spheres; i++)
  {
      Sphere sphere = spheres[i];

      Intersection intersectionSphere = CalculateSphereIntersection(shadowRay, sphere);
      if(intersectionSphere.object != nullptr && lightPosition.dist(intersection.position) > intersectionSphere.position.dist(intersection.position)) return true;
  }
  
  
  for(int i = 0; i < num_triangles; i++)
  {
      Triangle triangle = triangles[i];

      Intersection intersectionTriangle = CalculateTriangleIntersection(shadowRay, triangle);

      if(intersectionTriangle.object != nullptr && lightPosition.dist(intersection.position) > intersectionTriangle.position.dist(intersection.position)) return true;    
  }
  
  return false;
}

void ColorPixel(unsigned int row, unsigned int col, const Intersection& intersection, bool shouldPrintStuff = false)
{
  float red = 0.0f;
  float green = 0.0f;
  float blue = 0.0f;

  // Adding point lights
  for(int i = 0; i < num_lights; i++)
  {
    Light light = lights[i];

    if(!CheckInShadow(intersection, light, false))
    {
      Vector3 lightPos = Vector3(light.position);

      Vector3 toLight = lightPos - intersection.position;
      toLight.normalize();

      float lightIntersectionDotProduct = toLight.dot(intersection.normal);
      if (lightIntersectionDotProduct < 0.0f) { lightIntersectionDotProduct = 0.0f; }

      Vector3 reflection = toLight - intersection.normal * 2.0f * lightIntersectionDotProduct;
      reflection.normalize();

      float reflectionIntersectionDotProduct = reflection.dot(intersection.toCamera);
      if (reflectionIntersectionDotProduct < 0.0f) { reflectionIntersectionDotProduct = 0.0f; }
      
      red += light.color[0] * (intersection.color_diffuse.x * lightIntersectionDotProduct + intersection.color_specular.x * pow(reflectionIntersectionDotProduct, intersection.shininess));
      green += light.color[1] * (intersection.color_diffuse.y * lightIntersectionDotProduct + intersection.color_specular.y * pow(reflectionIntersectionDotProduct, intersection.shininess));
      blue += light.color[2] * (intersection.color_diffuse.z * lightIntersectionDotProduct + intersection.color_specular.z * pow(reflectionIntersectionDotProduct, intersection.shininess));    
    }
  }

  // Adding ambient light
  red += ambient_light[0];
  green += ambient_light[1];
  blue += ambient_light[2];

  red = glm::clamp(red, 0.0f, 1.0f);
  green = glm::clamp(green, 0.0f, 1.0f);
  blue = glm::clamp(blue, 0.0f, 1.0f);

  buffer[row][col][0] = (unsigned char)(red * 255.0f);
  buffer[row][col][1] = (unsigned char)(green * 255.0f);
  buffer[row][col][2] = (unsigned char)(blue * 255.0f);
}

#pragma region Spheres

void CreateSphereIntersectionColors(const Ray& ray, Sphere& sphere, Intersection& intersection)
{
  intersection.toCamera = intersection.position - ray.origin;
  intersection.toCamera.normalize();

  Vector3 spherePosition = Vector3(sphere.position);

  intersection.normal = intersection.position - spherePosition;
  intersection.normal.normalize();

  intersection.color_diffuse = Vector3(sphere.color_diffuse);
  intersection.color_specular = Vector3(sphere.color_specular);
  intersection.shininess = sphere.shininess;
}

Intersection CalculateSphereIntersection(const Ray& ray, Sphere& sphere)
{
  Vector3 sphereCenterToRayOrigin = ray.origin - Vector3(sphere.position);

  float b = 2.0f * ray.direction.dot(sphereCenterToRayOrigin);             
  float c = pow(sphereCenterToRayOrigin.x, 2.0) + pow(sphereCenterToRayOrigin.y, 2.0) + pow(sphereCenterToRayOrigin.z, 2.0) - pow(sphere.radius, 2.0);

  if(b * b <= 4.0f * c) return Intersection();

  float t0 = (-b - sqrt(b * b - 4.0f * c)) / 2.0f;
  float t1 = (-b + sqrt(b * b - 4.0f * c)) / 2.0f;
  float tMin = (t1 < t0) ? t1 : t0;

  if(tMin >= 0.0f)
  {
    Intersection intersection = Intersection();
    
    intersection.object = static_cast<Object*>(&sphere);
    intersection.position = ray.origin + ray.direction * tMin;

    return intersection;
  }
  
  return Intersection();
}

#pragma endregion

#pragma region Triangles

void CreateTriangleIntersectionColors(const Ray& ray, const Triangle& triangle, Intersection& intersection)
{
    intersection.toCamera = intersection.position - ray.origin;
    intersection.toCamera.normalize();

    Vector3 vertexA(triangle.v[0].position);
    Vector3 vertexB(triangle.v[1].position);
    Vector3 vertexC(triangle.v[2].position);

    Vector3 edgeA = vertexB - vertexA;
    Vector3 edgeB = vertexC - vertexA;
    Vector3 normal = edgeA.cross(edgeB);
    normal.normalize();

    Vector3 vAvB = vertexB - vertexA;
    Vector3 vBvC = vertexC - vertexB;
    Vector3 vCvA = vertexA - vertexC;

    Vector3 vAvC = vertexC - vertexA;

    Vector3 vAp = intersection.position - vertexA;
    Vector3 vBp = intersection.position - vertexB;
    Vector3 vCp = intersection.position - vertexC;

    if (normal.dot(vBvC.cross(vBp))  <= 0.0f) std::cout << "first if was: " << normal.dot(vBvC.cross(vBp)) << std::endl;
    if (normal.dot(vCvA.cross(vCp))  <= 0.0f) std::cout << "second if was: " << normal.dot(vCvA.cross(vCp))<< std::endl;
    if (normal.dot(vAvB.cross(vAp)) <= 0.0f) std::cout << "second if was: " << normal.dot(vAvB.cross(vAp)) << std::endl;
    
    float totalArea = (vAvB).cross(vAvC).length();
      
    // lambda1: Area of triangle BCP / totalArea
    float lambda1 = normal.dot(vBvC.cross(vBp)) / totalArea;

    // lambda2: Area of triangle CAP / totalArea
    float lambda2 = normal.dot(vCvA.cross(vCp)) / totalArea;

    // lambda3: Area of triangle ABP / totalArea
    float lambda3 = 1.0f - lambda1 - lambda2;
    
    Vector3 normalA(triangle.v[0].normal);
    Vector3 normalB(triangle.v[1].normal);
    Vector3 normalC(triangle.v[2].normal);
    intersection.normal = normalA * lambda1 + normalB * lambda2 + normalC * lambda3;
    intersection.normal.normalize();

    Vector3 diffuseA(triangle.v[0].color_diffuse);
    Vector3 diffuseB(triangle.v[1].color_diffuse);
    Vector3 diffuseC(triangle.v[2].color_diffuse);
    intersection.color_diffuse = diffuseA * lambda1 + diffuseB * lambda2 + diffuseC * lambda3;

    Vector3 specularA(triangle.v[0].color_specular);
    Vector3 specularB(triangle.v[1].color_specular);
    Vector3 specularC(triangle.v[2].color_specular);
    intersection.color_specular = specularA * lambda1 + specularB * lambda2 + specularC * lambda3;

    float shininessA = triangle.v[0].shininess;
    float shininessB = triangle.v[1].shininess;
    float shininessC = triangle.v[2].shininess;
    intersection.shininess = shininessA * lambda1 + shininessB * lambda2 + shininessC * lambda3;  
}

Intersection CalculateTriangleIntersection(const Ray &ray, Triangle &triangle)
{
    Vector3 vertexA(triangle.v[0].position);
    Vector3 vertexB(triangle.v[1].position);
    Vector3 vertexC(triangle.v[2].position);

    Vector3 edgeA = vertexB - vertexA;
    Vector3 edgeB = vertexC - vertexA;
    Vector3 normal = edgeA.cross(edgeB);
    normal.normalize();

    float denominator = normal.dot(ray.direction);
    if (fabs(denominator) < 1e-8) return Intersection();

    float d = -normal.dot(vertexA);
    float t = -(normal.dot(ray.origin) + d) / denominator;

    if (t <= 1e-8) return Intersection();

    Vector3 p = ray.origin + ray.direction * t;

    // inside-outside test
    Vector3 vAvB = vertexB - vertexA;
    Vector3 vBvC = vertexC - vertexB;
    Vector3 vCvA = vertexA - vertexC;

    Vector3 vAp = p - vertexA;
    Vector3 vBp = p - vertexB;
    Vector3 vCp = p - vertexC;

    if (normal.dot(vAvB.cross(vAp)) <= 1e-8) return Intersection();
    if (normal.dot(vBvC.cross(vBp)) <= 1e-8) return Intersection();
    if (normal.dot(vCvA.cross(vCp)) <= 1e-8) return Intersection();

    Intersection intersection = Intersection();

    intersection.object = static_cast<Object*>(&triangle);
    intersection.position = p;

    return intersection;
}

#pragma endregion

void display()
{
  float imageAspectRatio = WIDTH / (float)HEIGHT;
  Vector3 rayOrigin = Vector3(0.0f, 0.0f, 0.0f);

  for (unsigned int row = 0; row < HEIGHT; row++)
  {
     for (unsigned int col = 0; col < WIDTH; col++)
     {
        std::vector<Intersection> intersections;

        Vector3 origin = Vector3(0.0f, 0.0f, 0.0f);
        
        float Px = (2.0f * ((col + 0.5f) / WIDTH) - 1.0f) * tan(fov / 2.0 * M_PI / 180.0) * imageAspectRatio;
        float Py = (2.0f * ((row + 0.5f) / HEIGHT) - 1.0f) * tan(fov / 2.0 * M_PI / 180.0);
        Vector3 rayDirection = Vector3(Px, Py, -1.0f);

        Ray currentRay = Ray(rayOrigin, rayDirection);
        rays[row][col] = currentRay;

        for (int i = 0; i < num_spheres; i++)
        {
            Sphere sphere = spheres[i];
            Intersection intersection = CalculateSphereIntersection(currentRay, sphere);
            if(intersection.object != nullptr)
            {
              CreateSphereIntersectionColors(currentRay, sphere, intersection);
              intersections.push_back(intersection);
            }
        }

        for (int i = 0; i < num_triangles; i++)
        {
            Triangle triangle = triangles[i];
            Intersection intersection = CalculateTriangleIntersection(currentRay, triangle);
            if(intersection.object != nullptr)
            {
              CreateTriangleIntersectionColors(currentRay, triangle, intersection);
              intersections.push_back(intersection);
            }
        }

        if (!intersections.empty())
        {
          // Find the closest intersection

          Intersection* closestIntersection = nullptr;
          float minDist = std::numeric_limits<float>::max();

          for (unsigned int i = 0; i < intersections.size(); i++) 
          {
            Intersection* intersection = &intersections[i];
            if (intersection != nullptr)
            {
              float distFromCamera = rayOrigin.dist(intersection->position);

              if (distFromCamera < minDist)
              {
                  minDist = distFromCamera;
                  closestIntersection = intersection;
              }
            }
          }

          // Color the closest intersection and clear intersections
          if (closestIntersection != nullptr)
          {
            ColorPixel(row, col, *closestIntersection);
          }
        }
        else // Draw Background if no intersections exist
        {
          buffer[row][col][0] = (unsigned char)(255.0f);
          buffer[row][col][1] = (unsigned char)(255.0f);
          buffer[row][col][2] = (unsigned char)(255.0f);
        }
     }
  }
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // Hack to make it only draw once.
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

