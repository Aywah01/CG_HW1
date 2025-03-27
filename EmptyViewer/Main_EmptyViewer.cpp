#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>
#include <memory>
#include <limits>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// Global variables for image dimensions
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;

// Ray class represents a ray with origin and direction
class Ray {
public:
    vec3 origin;
    vec3 direction;

    Ray(const vec3& o, const vec3& d) : origin(o), direction(d) {}
};

// -------------------------------------------------
// Surface Class (Abstract Base Class)
// -------------------------------------------------
class Surface {
public:
    virtual bool intersect(const Ray& ray, float& t) const = 0;
};

// -------------------------------------------------
// Plane Class
// -------------------------------------------------
class Plane : public Surface {
public:
    float y;

    Plane(float y) : y(y) {}

    bool intersect(const Ray& ray, float& t) const override {
        // Check for parallel rays
        if (ray.direction.y == 0) return false; // Parallel to the plane
        t = (y - ray.origin.y) / ray.direction.y;
        return t >= 0; // Ensure intersection is in front of the camera
    }
};

// -------------------------------------------------
// Sphere Class
// -------------------------------------------------
class Sphere : public Surface {
public:
    vec3 center;
    float radius;

    Sphere(const vec3& c, float r) : center(c), radius(r) {}

    bool intersect(const Ray& ray, float& t) const override {
        vec3 oc = ray.origin - center;
        float a = dot(ray.direction, ray.direction);
        float b = 2.0f * dot(oc, ray.direction);
        float c = dot(oc, oc) - radius * radius;
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0) return false; // No intersection

        float sqrt_disc = glm::sqrt(discriminant);
        float t0 = (-b - sqrt_disc) / (2.0f * a);
        float t1 = (-b + sqrt_disc) / (2.0f * a);

        t = (t0 < t1) ? t0 : t1; // Choose closest positive intersection
        return t >= 0;
    }
};

// -------------------------------------------------
// Camera Class
// -------------------------------------------------
class Camera {
public:
    vec3 eye;
    vec3 u, v, w; // Orientation vectors
    float l, r, b, t, d; // View frustum parameters
    int nx, ny; // Image resolution

    Camera(const vec3& e, const vec3& u, const vec3& v, const vec3& w, float l, float r, float b, float t, float d, int nx, int ny)
        : eye(e), u(u), v(v), w(w), l(l), r(r), b(b), t(t), d(d), nx(nx), ny(ny) {
    }

    Ray generateRay(float i, float j) const {
        float u_coord = l + (r - l) * (i + 0.5f) / nx;
        float v_coord = b + (t - b) * (j + 0.5f) / ny;
        vec3 direction = normalize(-d * w + u_coord * u + v_coord * v);
        return Ray(eye, direction);
    }
};

// -------------------------------------------------
// Scene Class
// -------------------------------------------------
class Scene {
public:
    std::vector<Surface*> objects;

    void add(Surface* object) {
        objects.push_back(object);
    }

    bool intersect(const Ray& ray, float& t) const {
        float t_min = std::numeric_limits<float>::max();
        bool hit = false;
        for (const auto& obj : objects) {
            float t_obj;
            if (obj->intersect(ray, t_obj) && t_obj < t_min) {
                t_min = t_obj;
                hit = true;
            }
        }
        t = t_min;
        return hit;
    }
};

// -------------------------------------------------
// Global Scene and Camera
// -------------------------------------------------
Scene scene;
Camera camera(vec3(0, 0, 0), vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1), -0.1f, 0.1f, -0.1f, 0.1f, 0.1f, 512, 512);

// -------------------------------------------------
// Render Function
// -------------------------------------------------
void render() {
    OutputImage.clear();
    for (int j = 0; j < Height; ++j) {
        for (int i = 0; i < Width; ++i) {
            Ray ray = camera.generateRay(i, j);
            float t;
            if (scene.intersect(ray, t)) {
                OutputImage.push_back(1.0f); // R
                OutputImage.push_back(1.0f); // G
                OutputImage.push_back(1.0f); // B
            }
            else {
                OutputImage.push_back(0.0f); // R
                OutputImage.push_back(0.0f); // G
                OutputImage.push_back(0.0f); // B
            }
        }
    }
}


// -------------------------------------------------
// Resize Callback
// -------------------------------------------------
void resize_callback(GLFWwindow*, int nw, int nh) {
    Width = nw;
    Height = nh;
    glViewport(0, 0, nw, nh);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, static_cast<double>(Width), 0.0, static_cast<double>(Height), 1.0, -1.0);
    OutputImage.reserve(Width * Height * 3);
    render();
}

// -------------------------------------------------
// Main Function
// -------------------------------------------------
int main(int argc, char* argv[]) {
    // Initialize Scene and Add objects to the scene
    scene.add(new Plane(-2.0f));
    scene.add(new Sphere(vec3(-4, 0, -7), 1.0f));
    scene.add(new Sphere(vec3(0, 0, -7), 2.0f));
    scene.add(new Sphere(vec3(4, 0, -7), 1.0f));

    // Initialize Window
    GLFWwindow* window;

    if (!glfwInit())
        return -1;

    window = glfwCreateWindow(Width, Height, "Simple Ray Tracer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glfwSetFramebufferSizeCallback(window, resize_callback);
    resize_callback(NULL, Width, Height);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS || glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}