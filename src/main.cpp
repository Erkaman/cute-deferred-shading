#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cstdlib>
#include <cstring>
#include <string>
#include <chrono>
#include <ctime>
#include <thread>
#include <algorithm>

#include "lodepng.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using std::string;
using std::vector;

//
//
// Begin OpenGL utility.
//
//

inline void CheckOpenGLError(const char* stmt, const char* fname, int line)
{
    GLenum err = glGetError();
    if (err != GL_NO_ERROR){
        printf("OpenGL error %08x, at %s:%i - for %s.\n", err, fname, line, stmt);
        exit(1);
    }
}

// helper macro that checks for GL errors.
#define GL_C(stmt) do {					\
	stmt;						\
	CheckOpenGLError(#stmt, __FILE__, __LINE__);	\
    } while (0)

inline char* GetShaderLogInfo(GLuint shader) {
    GLint len;
    GLsizei actualLen;
    GL_C(glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len));
    char* infoLog = new char[len];
    GL_C(glGetShaderInfoLog(shader, len, &actualLen, infoLog));
    return infoLog;
}

inline GLuint CreateShaderFromString(const std::string& shaderSource, const GLenum shaderType) {
    GLuint shader;

    GL_C(shader = glCreateShader(shaderType));
    const char *c_str = shaderSource.c_str();
    GL_C(glShaderSource(shader, 1, &c_str, NULL));
    GL_C(glCompileShader(shader));

    GLint compileStatus;
    GL_C(glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus));
    if (compileStatus != GL_TRUE) {
        printf("Could not compile shader\n\n%s \n\n%s\n", shaderSource.c_str(),
               GetShaderLogInfo(shader));
        exit(1);
    }

    return shader;
}

/*
  Load shader with only vertex and fragment shader.
*/
inline GLuint LoadNormalShader(const std::string& vsSource, const std::string& fsShader){
    GLuint vs = CreateShaderFromString(vsSource, GL_VERTEX_SHADER);
    GLuint fs = CreateShaderFromString(fsShader, GL_FRAGMENT_SHADER);

    GLuint shader = glCreateProgram();
    glAttachShader(shader, vs);
    glAttachShader(shader, fs);
    glLinkProgram(shader);

    GLint Result;
    glGetProgramiv(shader, GL_LINK_STATUS, &Result);
    if (Result == GL_FALSE) {
        printf("Could not link shader \n\n%s\n", GetShaderLogInfo(shader));
        exit(1);
    }

    glDetachShader(shader, vs);
    glDetachShader(shader, fs);

    glDeleteShader(vs);
    glDeleteShader(fs);

    return shader;
}


//
//
// Begin vec3.
//
//

class vec3 {
public:
    float x, y, z;

    vec3(float x, float y, float z) { this->x = x; this->y = y; this->z = z; }

    vec3() { this->x = this->y = this->z = 0; }

    vec3& operator+=(const vec3& b) { (*this) = (*this) + b; return (*this); }

    friend vec3 operator-(const vec3& a, const vec3& b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
    friend vec3 operator+(const vec3& a, const vec3& b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
    friend vec3 operator*(const float s, const vec3& a) { return vec3(s * a.x, s * a.y, s * a.z); }
    friend vec3 operator*(const vec3& a, const float s) { return s * a; }

    static float length(const vec3& a) { return sqrt(vec3::dot(a, a)); }

    // dot product.
    static float dot(const vec3& a, const vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

    static float distance(const vec3& a, const vec3& b) { return length(a - b); }
    static vec3 normalize(const vec3& a) { return (1.0f / vec3::length(a)) * a; }

    // cross product.
    static vec3 cross(const vec3& a, const vec3& b) { return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

    //
    // Rotate the vector 'v' around the 'axis' for 'theta' degrees.
    // This is basically Rodrigues' rotation formula.
    //
    static vec3 rotate(const vec3& v, const float theta, const vec3& axis) {
        vec3 k = vec3::normalize(axis); // normalize for good measure.
        return v * cos(theta) + vec3::cross(k, v)* sin(theta) + (k * vec3::dot(k, v)) * (1.0f - cos(theta));
    }
};

//
//
// Begin mat4
//
//

class mat4 {
public:
    float m[4][4];

    mat4() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m[i][j] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    // return perspective projection matrix.
    static mat4 perspective(float fovy, float aspect, float zNear, float zFar) {
        const float tanHalfFovy = tan(fovy / 2.0f);
        mat4 m;
        m.m[0][0] = 1 / (aspect * tanHalfFovy);
        m.m[1][1] = 1 / (tanHalfFovy);
        m.m[2][2] = -(zFar + zNear) / (zFar - zNear);
        m.m[2][3] = -1;
        m.m[3][2] = -(2 * zFar * zNear) / (zFar - zNear);
        return m;
    }

    static mat4 lookAt(const vec3& eye, const vec3& center, const vec3& up)
        {
            mat4 m;
            // compute the basis vectors.
            vec3 forward = vec3::normalize(center - eye); // forward vector.
            vec3 left = vec3::normalize(vec3::cross(forward, up)); // left vector.
            vec3 u = vec3::cross(left, forward); // up vector.

            m.m[0][0] = left.x;
            m.m[1][0] = left.y;
            m.m[2][0] = left.z;
            m.m[0][1] = u.x;
            m.m[1][1] = u.y;
            m.m[2][1] = u.z;
            m.m[0][2] = -forward.x;
            m.m[1][2] = -forward.y;
            m.m[2][2] = -forward.z;
            m.m[3][0] = -vec3::dot(left, eye);
            m.m[3][1] = -vec3::dot(u, eye);
            m.m[3][2] = vec3::dot(forward, eye);
            return m;
        }

    friend mat4 operator*(const mat4& a, const mat4& b) {
        mat4 m;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m.m[i][j] = 0.0f;
                for (int k = 0; k < 4; k++) {
                    m.m[i][j] += a.m[i][k] * b.m[k][j];
                }
            }
        }
        return m;
    }
};


vec3 lerp(const vec3& a, const vec3& b, float t) {
    return vec3(
        a.x * (1.0 - t) + b.x * t,
        a.y * (1.0 - t) + b.y * t,
        a.z * (1.0 - t) + b.z * t

        );
}

//
//
// begin Camera
//
//
class Camera{
private:

    vec3 viewDir;
    vec3 right;
    vec3 up;
    vec3 position;

    double prevMouseX = 0.0;
    double prevMouseY = 0.0;

    double curMouseX = 0.0;
    double curMouseY = 0.0;

public:
    mat4 GetViewMatrix() {
        return mat4::lookAt(position, position + viewDir, up);

    }

    Camera(const vec3& position_, const vec3& viewDir_) : position(position_), viewDir(viewDir_){
        viewDir = vec3::normalize(viewDir);
        right = vec3::normalize(vec3::cross(viewDir, vec3(0.0f, 1.0f, 0.0f)));
        up = vec3(0, 1.0f, 0);
    }

    void Update(const float delta, GLFWwindow* window) {
        // we use mouse movement to change the camera viewing angle.
        prevMouseX = curMouseX;
        prevMouseY = curMouseY;
        glfwGetCursorPos(window, &curMouseX, &curMouseY);
        float mouseDeltaX = (float)(curMouseX - prevMouseX);
        float mouseDeltaY = (float)(curMouseY - prevMouseY);

        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            viewDir = vec3::normalize(vec3::rotate(viewDir, mouseDeltaX*-0.01f, up));
            viewDir = vec3::normalize(vec3::rotate(viewDir, mouseDeltaY*-0.01f, right));
            right = vec3::normalize(vec3::cross(viewDir, vec3(0.0f, 1.0f, 0.0f)));
            up = vec3(0.0f, 1.0f, 0.0f);
        }

        static float cameraSpeed;

        if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
            cameraSpeed = 1000.0f;
        else
            cameraSpeed = 100.0f;

        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            position += delta * cameraSpeed * viewDir;
        else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            position += -delta * cameraSpeed * viewDir;

        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            position += -delta * cameraSpeed * right;
        else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            position += +delta * cameraSpeed * right;

        if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
            position += +delta * cameraSpeed * up;
        else if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
            position += -delta * cameraSpeed * up;

        if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS) {
            printf("vec3(%f,%f,%f), vec3(%f,%f,%f), vec3(%f,%f,%f), vec3(%f,%f,%f),\n\n",
                   viewDir.x, viewDir.y, viewDir.z,
                   right.x, right.y, right.z,
                   up.x, up.y, up.z,
                   position.x, position.y, position.z
                );
        }

    }

    vec3 GetPosition() const {
        return position;
    }
};

//
//
// Now the actual demo starts.
//
//
// first comes global variables:
struct Mesh;
struct Material;
std::vector<Mesh*> meshes; // all the meshes of the .obj model.
std::vector<Material> materials; // all the materials of the .obj model.

// these four textures are the gbuffer.
GLuint colorTexture;
GLuint depthRenderbuffer;
GLuint normalTexture;
GLuint positionTexture;
// we use this fbo to render to the gbuffer.
GLuint fbo;

GLuint vao;

const int WINDOW_WIDTH = 1497;
const int WINDOW_HEIGHT = 1014;
GLFWwindow* window;
int FRAME_RATE = 60;
float totalTime = 0.0f; // global time.

int fbWidth, fbHeight;

GLuint outputGeoShader;
GLuint directionalLightShader;
GLuint pointLightShader;

// light sphere geometry:
GLuint spherePositionVbo;
GLuint sphereIndexVbo;
GLuint sphereIndexCount;

std::string sponzaDir; // directory where sponza is located.


Camera camera(vec3(0.0f, 0.0f, 0.0f), vec3::normalize(vec3(0.3f, 0.5f, 0.3f)));

struct Material {
    // to keep things simple, we ignore everything but the diffuse texture.
    GLuint diffuseTex;
    std::string diffuseTexFile;
    std::string name;
};

struct Mesh {
    GLuint positionVbo;
    GLuint normalVbo;
    GLuint uvVbo;

    std::string name;

    int matId;

    std::vector<float> positions;
    std::vector<float> normals;
    std::vector<float> uvs;

    // upload mesh to GPU.
    void UploadMesh() {
        GL_C(glGenBuffers(1, &this->positionVbo));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, this->positionVbo));
        GL_C(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*this->positions.size(), this->positions.data(), GL_STATIC_DRAW));

        GL_C(glGenBuffers(1, &this->normalVbo));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, this->normalVbo));
        GL_C(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*this->normals.size(), this->normals.data(), GL_STATIC_DRAW));

        GL_C(glGenBuffers(1, &this->uvVbo));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, this->uvVbo));
        GL_C(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*this->uvs.size(), this->uvs.data(), GL_STATIC_DRAW));
    }
};

GLuint LoadTexture(const char* file) {
    std::vector<unsigned char> buffer;
    lodepng::load_file(buffer, file);

    lodepng::State state;
    unsigned int width, height;
    std::vector<unsigned char> imageData;
    unsigned error = lodepng::decode(imageData, width, height, state, buffer);

    if (error != 0) {
        printf("Could not load texture %s: %s\n", file, lodepng_error_text(error));
    }

    GLuint tex;
    GL_C(glGenTextures(1, &tex));

    GL_C(glBindTexture(GL_TEXTURE_2D, tex));

    GL_C(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, imageData.data()));

    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT));
    GL_C(glGenerateMipmap(GL_TEXTURE_2D));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR));

    GL_C(glBindTexture(GL_TEXTURE_2D, 0));

    return tex;
}

// replace '\' with '/' in filename, so that we consistently use '/'
string FixFilename(string file) {
    std::replace(file.begin(), file.end(), '\\', '/');
    return file;
}

std::string strReplace(std::string &s,
                       const std::string &toReplace,
                       const std::string &replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

void LoadModel(void) {
    sponzaDir = FixFilename(sponzaDir);
    if (sponzaDir[sponzaDir.size() - 1] != '/') {
        sponzaDir += "/";
    }

    std::string mtlDir = sponzaDir;
    std::string inputfile = sponzaDir;
    inputfile += "sponza.obj";

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> objMaterials;
    std::string err;
    bool ret = tinyobj::LoadObj(&attrib, &shapes, &objMaterials, &err, inputfile.c_str(), mtlDir.c_str(), true); // do triangulate.

    if (!err.empty()) {
        printf("FAILED LOAD OBJ MODEL: %s\n", err.c_str());
    }
    if (!ret) {
        exit(1);
    }

    // get material info.
    for (size_t i = 0; i < objMaterials.size(); i++) {
        tinyobj::material_t& m = objMaterials[i];

        Material newMat;
        newMat.name = m.name;

        if (m.diffuse_texname == "") {
            newMat.diffuseTexFile = "";
        }
        else {
            newMat.diffuseTexFile = FixFilename(mtlDir + m.diffuse_texname);
        }
        materials.push_back(newMat);
    }

    for (size_t i = 0; i < shapes.size(); i++) {
        tinyobj::shape_t& s = shapes[i];

        tinyobj::mesh_t m = s.mesh;

        Mesh* mesh = new Mesh();
        mesh->matId = m.material_ids[0];

        for (size_t j = 0; j < m.indices.size(); j += 3) {
            if (mesh->matId != m.material_ids[j / 3]) {
                // new material. so start a new mesh.
                mesh->UploadMesh();
                meshes.push_back(mesh);

                mesh = new Mesh();
                mesh->matId = m.material_ids[j / 3];
            }

            tinyobj::index_t i0 = m.indices[j + 0];
            tinyobj::index_t i1 = m.indices[j + 1];
            tinyobj::index_t i2 = m.indices[j + 2];

            // pos
            mesh->positions.push_back(attrib.vertices[i0.vertex_index * 3 + 0]);
            mesh->positions.push_back(attrib.vertices[i0.vertex_index * 3 + 1]);
            mesh->positions.push_back(attrib.vertices[i0.vertex_index * 3 + 2]);

            // normal.
            mesh->normals.push_back(attrib.normals[i0.normal_index * 3 + 0]);
            mesh->normals.push_back(attrib.normals[i0.normal_index * 3 + 1]);
            mesh->normals.push_back(attrib.normals[i0.normal_index * 3 + 2]);

            // uvs
            mesh->uvs.push_back(attrib.texcoords[i0.texcoord_index * 2 + 0]);
            mesh->uvs.push_back(attrib.texcoords[i0.texcoord_index * 2 + 1]);

            // pos
            mesh->positions.push_back(attrib.vertices[i1.vertex_index * 3 + 0]);
            mesh->positions.push_back(attrib.vertices[i1.vertex_index * 3 + 1]);
            mesh->positions.push_back(attrib.vertices[i1.vertex_index * 3 + 2]);

            // normal.
            mesh->normals.push_back(attrib.normals[i1.normal_index * 3 + 0]);
            mesh->normals.push_back(attrib.normals[i1.normal_index * 3 + 1]);
            mesh->normals.push_back(attrib.normals[i1.normal_index * 3 + 2]);

            // uvs
            mesh->uvs.push_back(attrib.texcoords[i1.texcoord_index * 2 + 0]);
            mesh->uvs.push_back(attrib.texcoords[i1.texcoord_index * 2 + 1]);

            // pos
            mesh->positions.push_back(attrib.vertices[i2.vertex_index * 3 + 0]);
            mesh->positions.push_back(attrib.vertices[i2.vertex_index * 3 + 1]);
            mesh->positions.push_back(attrib.vertices[i2.vertex_index * 3 + 2]);

            // normal.
            mesh->normals.push_back(attrib.normals[i2.normal_index * 3 + 0]);
            mesh->normals.push_back(attrib.normals[i2.normal_index * 3 + 1]);
            mesh->normals.push_back(attrib.normals[i2.normal_index * 3 + 2]);

            // uvs
            mesh->uvs.push_back(attrib.texcoords[i2.texcoord_index * 2 + 0]);
            mesh->uvs.push_back(attrib.texcoords[i2.texcoord_index * 2 + 1]);
        }

        mesh->name = s.name;

        meshes.push_back(mesh);
    }

    for (Material& m : materials) {
        if (m.diffuseTexFile != "") {
            m.diffuseTex = LoadTexture(m.diffuseTexFile.c_str());
        }
    }

    for (Mesh* mesh : meshes) {
        mesh->UploadMesh();
    }
}

void InitGlfw() {
    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Deferred Shading Demo", NULL, NULL);
    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(window);

    // load GLAD.
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    // Bind and create VAO, otherwise, we can't do anything in OpenGL.
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);


    glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
}

// create simple UV-sphere.
void CreateSphere() {
    int stacks = 20;
    int slices = 20;
    const float PI = 3.14f;

    std::vector<float> positions;
    std::vector<GLuint> indices;

    // loop through stacks.
    for (int i = 0; i <= stacks; ++i){

        float V = (float)i / (float)stacks;
        float phi = V * PI;

        // loop through the slices.
        for (int j = 0; j <= slices; ++j){

            float U = (float)j / (float)slices;
            float theta = U * (PI * 2);

            // use spherical coordinates to calculate the positions.
            float x = cos(theta) * sin(phi);
            float y = cos(phi);
            float z = sin(theta) * sin(phi);

            positions.push_back(x);
            positions.push_back(y);
            positions.push_back(z);
        }
    }

    // Calc The Index Positions
    for (int i = 0; i < slices * stacks + slices; ++i){
        indices.push_back(GLuint(i));
        indices.push_back(GLuint(i + slices + 1));
        indices.push_back(GLuint(i + slices));

        indices.push_back(GLuint(i + slices + 1));
        indices.push_back(GLuint(i));
        indices.push_back(GLuint(i + 1));
    }

    // upload geometry to GPU.
    GL_C(glGenBuffers(1, &spherePositionVbo));
    GL_C(glBindBuffer(GL_ARRAY_BUFFER, spherePositionVbo));
    GL_C(glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*positions.size(), positions.data(), GL_STATIC_DRAW));

    GL_C(glGenBuffers(1, &sphereIndexVbo));
    GL_C(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereIndexVbo));
    GL_C(glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*indices.size(), indices.data(), GL_STATIC_DRAW));

    sphereIndexCount = indices.size();
}

// configure a shader for usage in deferred rendering.
void SetupDeferredShader(GLuint shader) {
    // bind gbuffer textures.

    GL_C(glActiveTexture(GL_TEXTURE0 + 0));
    GL_C(glBindTexture(GL_TEXTURE_2D, colorTexture));
    GL_C(glUniform1i(glGetUniformLocation(shader, "uColorTex"), 0));

    GL_C(glActiveTexture(GL_TEXTURE0 + 1));
    GL_C(glBindTexture(GL_TEXTURE_2D, normalTexture));
    GL_C(glUniform1i(glGetUniformLocation(shader, "uNormalTex"), 1));

    GL_C(glActiveTexture(GL_TEXTURE0 + 2));
    GL_C(glBindTexture(GL_TEXTURE_2D, positionTexture));
    GL_C(glUniform1i(glGetUniformLocation(shader, "uPositionTex"), 2));

    GL_C(glUniform3f(glGetUniformLocation(shader, "uCameraPos"), camera.GetPosition().x, camera.GetPosition().y, camera.GetPosition().z));
}

void RenderPointLight(float radius, const vec3& position, const vec3& color) {
    GL_C(glUniform1f(glGetUniformLocation(pointLightShader, "uLightRadius"), radius));
    GL_C(glUniform3f(glGetUniformLocation(pointLightShader, "uLightPosition"), position.x, position.y, position.z));
    GL_C(glUniform3f(glGetUniformLocation(pointLightShader, "uLightColor"), color.x, color.y, color.z));
    GL_C(glDrawElements(GL_TRIANGLES, sphereIndexCount, GL_UNSIGNED_INT, 0));
}

void Render() {
    // setup matrices.
    mat4 projectionMatrix = mat4::perspective(0.9f, (float)(WINDOW_WIDTH) / WINDOW_HEIGHT, 0.1f, 3000.0f);
    mat4 viewMatrix = camera.GetViewMatrix();
    mat4 VP = viewMatrix * projectionMatrix;

    // setup GL state.
    GL_C(glEnable(GL_DEPTH_TEST));
    GL_C(glDepthMask(true));
    GL_C(glDisable(GL_BLEND));
    GL_C(glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE));
    GL_C(glEnable(GL_CULL_FACE));
    GL_C(glFrontFace(GL_CCW));

    //
    // In the first pass, we just write to the gbuffer.
    //
    GL_C(glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo)); // bind g buffer for writing.

    GL_C(glViewport(0, 0, fbWidth, fbHeight));
    GL_C(glClearColor(0.0f, 0.0f, 0.3f, 1.0f));
    GL_C(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));

    GL_C(glUseProgram(outputGeoShader)); // "output geometry to gbuffer" shader
    GL_C(glUniformMatrix4fv(glGetUniformLocation(outputGeoShader, "uVp"), 1, GL_FALSE, (GLfloat *)VP.m));

    // now we render all the meshes, one after one.
    for (Mesh* mesh : meshes) {
        //
        // setup vertex attribs.
        //
        GL_C(glEnableVertexAttribArray(0));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, mesh->positionVbo));
        GL_C(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0));

        GL_C(glEnableVertexAttribArray(1));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, mesh->normalVbo));
        GL_C(glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0));

        GL_C(glEnableVertexAttribArray(2));
        GL_C(glBindBuffer(GL_ARRAY_BUFFER, mesh->uvVbo));
        GL_C(glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, (void*)0));

        // enable texture.

        Material mat = materials[mesh->matId];
        if (mat.diffuseTexFile != "") {
            GL_C(glActiveTexture(GL_TEXTURE0));
            GL_C(glBindTexture(GL_TEXTURE_2D, mat.diffuseTex));
	    GL_C(glUniform1i(glGetUniformLocation(outputGeoShader, "uDiffTex"), 0));
        }

        glDrawArrays(GL_TRIANGLES, 0, mesh->positions.size());
    }

    GL_C(glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0)); // stop writing to gbuffer.

    //
    // Now comes the Deferred shading!
    //
    GL_C(glViewport(0, 0, fbWidth, fbHeight));
    GL_C(glClearColor(0.0f, 0.0f, 0.3f, 1.0f));
    GL_C(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));

    //
    // first, we render a single directional light, with a fullscreen pass.
    //
    GL_C(glUseProgram(directionalLightShader));
    SetupDeferredShader(directionalLightShader);
    // we use attribute-less rendering to render a full-screen triangle.
    // so the triangle vertices are basically stored in the vertex shader.
    // see the vertex shader for more details.
    GL_C(glDrawArrays(GL_TRIANGLES, 0, 3));

    //
    // Next, we render all the point light soures.
    // We will be doing our own depth testing in frag shader, so disable depth testing.
    // Enable alpha blending. So that the rendered point lights are added to the framebuffer.
    //
    GL_C(glDisable(GL_DEPTH_TEST));
    GL_C(glEnable(GL_BLEND));
    GL_C(glBlendFunc(GL_ONE, GL_ONE));

    // We render only the inner faces of the light sphere.
    // In other words, we render the back-faces and not the front-faces of the sphere.
    // If we render the front-faces, the lighting of the light sphere disappears if
    // we are inside the sphere, which is weird. But by rendering the back-faces instead,
    // we solve this problem.
    GL_C(glFrontFace(GL_CW));

    GL_C(glUseProgram(pointLightShader));
    SetupDeferredShader(pointLightShader);
    GL_C(glUniformMatrix4fv(glGetUniformLocation(pointLightShader, "uVp"), 1, GL_FALSE, (GLfloat*)VP.m));
    // We render every point light as a light sphere. And this light sphere is added onto the framebuffer
    // with additive alpha blending.
    GL_C(glEnableVertexAttribArray(0));
    GL_C(glBindBuffer(GL_ARRAY_BUFFER, spherePositionVbo));
    GL_C(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0));
    GL_C(glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereIndexVbo));


    //
    // Now it's time to render a TON of point lights.
    //

    const int PS = 14;
    vec3 palette[PS] = {
        vec3(1.0f, 0.0f, 0.0f),
        vec3(0.0f, 1.0f, 0.0f),
        vec3(0.0f, 0.0f, 1.0f),
        vec3(1.0f, 1.0f, 0.0f),
        vec3(1.0f, 0.0f, 1.0f),
        vec3(0.0f, 1.0f, 1.0f),
        vec3(1.0f, 1.0f, 1.0f),
        vec3(0.9f, 0.5f, +.3f),
        vec3(0.5f, 0.5f, 0.5f),
        vec3(0.3f, 0.6f, 0.9f),
        vec3(0.4f, 0.8f, 0.4f),
        vec3(0.3f, 0.3f, 1.0f),
        vec3(1.0f, 0.5f, 0.8f),
        vec3(0.5f, 0.9f, 0.2f),
    };
    const float PI = 3.14;

    int N;
    float RX;
    float RZ;

    float v;
    float n;
    float d;
    float k;

    N = 20;
    RX = 1000.0f;
    RZ = 800.0f;

    v = 0.6;
    n = 7;
    d = 4;
    k = n / d;
    // place the point lights on a rose curve:
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime*0.2;
        vec3 p = vec3(
            RX * cos(k * theta * v) * cos(theta * v),
            90,
            RZ * cos(k * theta * v) * sin(theta * v)
            );
        RenderPointLight(270.0f, p, palette[(i + 2) % PS]);
    }

    N = 20;
    RX = 900.0f;
    RZ = 400.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            490,
            RZ * sin(theta)
            );
        RenderPointLight(270.0f, p, palette[(i + 3) % PS]);
    }

    N = 20;
    RX = 900.0f;
    RZ = 300.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            250,
            RZ * sin(theta)
            );
        RenderPointLight(270.0f, p, palette[(i + 7) % PS]);
    }



    N = 20;
    RX = 1000.0f;
    RZ = 700.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            170,
            RZ * sin(theta)
            );
        RenderPointLight(350.0f, p, palette[(i + 0) % PS]);
    }


    N = 10;
    RX = 1000.0f;
    RZ = 700.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            290,
            RZ * sin(theta)
            );
        RenderPointLight(350.0f, p, palette[(i + 7) % PS]);
    }


    N = 20;
    RX = 900.0f;
    RZ = 400.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            690,
            RZ * sin(theta)
            );
        RenderPointLight(270.0f, p, palette[(i + 9) % PS]);
    }


    N = 30;
    RX = 1300.0f;
    RZ = 170.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            250,
            RZ * sin(theta)
            );
        RenderPointLight(270.0f, p, palette[(i + 7) % PS]);
    }

    N = 20;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            1250,
            600 * cos(theta),
            300 * sin(theta)
            );
        RenderPointLight(350, p, palette[(i + 7) % PS]);
    }


    N = 20;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            -1200,
            300 + 300 * cos(theta),
            300 * sin(theta)
            );
        RenderPointLight(350, p, palette[(i + 7) % PS]);
    }

    N = 20;
    RX = 800.0f;
    RZ = 140.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta) - 200,
            390,
            RZ * sin(theta)
            );
        RenderPointLight(350.0f, p, palette[(i + 10) % PS]);
    }

    N = 20;
    RX = 800.0f;
    RZ = 140.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta) - 150,
            890,
            RZ * sin(theta)
            );
        RenderPointLight(350.0f, p, palette[(i + 8) % PS]);
    }

    N = 20;
    RX = 800.0f;
    RZ = 140.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            1300,
            RZ * sin(theta)
            );
        RenderPointLight(350.0f, p, palette[(i + 2) % PS]);
    }

    N = 20;
    RX = 900.0f;
    RZ = 400.0f;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            RX * cos(theta),
            1300 + 100, //200.0f*sin(totalTime*0.3),
            RZ * sin(theta)
            );
        RenderPointLight(400.0f, p, palette[(i + 4) % PS]);
    }


    N = 20;
    for (int i = 0; i < N; i++) {
        float theta = 2.0f * (i / (float)N) * 2 * PI;
        theta += totalTime * 0.1f;
        vec3 p = vec3(
            //   -1400 + 100 * cos(totalTime * 0.2),
            -1421,
            340 + 300 * cos(theta),
            300 * sin(theta)
            );
        RenderPointLight(320, p, palette[(i + 7) % PS]);
    }

}

void HandleInput() {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    camera.Update(1.0 / (float)FRAME_RATE, window);
}

int main(int argc, char** argv) {
    if (argc == 1) {
        printf("Please give the directory of Sponza as a command line argument\n");
        exit(1);
    }
    else {
        sponzaDir = argv[1];
    }

    InitGlfw();

    outputGeoShader = LoadNormalShader(
        "#version 330\n"
        "layout(location = 0) in vec3 vsPos;"
        "layout(location = 1) in vec3 vsNormal;"
        "layout(location = 2) in vec2 vsUv;"

        "out vec3 fsPos;"
        "out vec3 fsNormal;"
        "out vec2 fsUv;"

        "uniform mat4 uVp;"
        "void main()"
        "{"
        "  fsPos = vsPos;"
        "  fsNormal = vsNormal;"
        "  fsUv = vsUv;"
        "  gl_Position = uVp * vec4(vsPos, 1.0);"
        "}",

        "#version 330\n"
        "in vec3 fsPos;"
        "in vec3 fsNormal;"
        "in vec2 fsUv;"

        "uniform sampler2D uDiffTex;"

        "out vec4 geoData[3];"

        "void main()"
        "{"
        // seems like the textures in Sponza are flipped. So flip then.
        "  vec4 diff = texture(uDiffTex, vec2(1.0, -1.0)*fsUv).rgba;"
        "  if (diff.a < 0.2) {"
        "    discard;"
        "  }"

        // output geometry.
        "  geoData[0] = vec4(diff.rgb, 1);"
        "  geoData[1] = vec4(fsNormal, 1);"
        "  geoData[2] = vec4(fsPos, 1);"
        "}"
        );


    directionalLightShader = LoadNormalShader(
        "#version 330\n"
        "out vec2 fsUv;"

        // full screen triangle vertices.
        "const vec2 verts[3] = vec2[](vec2(-1, -1), vec2(3, -1), vec2(-1, 3));"
        "const vec2 uvs[3] = vec2[](vec2(0, 0), vec2(2, 0), vec2(0, 2));"

        "void main()"
        "{"
        "  fsUv = uvs[gl_VertexID];"
        "  gl_Position = vec4( verts[gl_VertexID], 0.0, 1.0);"
        "}"
        ,

        "#version 330\n"
        "in vec2 fsUv;"

        "out vec4 outColor;"

        "uniform sampler2D uColorTex;"
        "uniform sampler2D uNormalTex;"
        "uniform sampler2D uPositionTex;"

        "uniform vec3 uCameraPos;"

        "void main()"
        "{"
        "  vec3 albedo = texture(uColorTex, fsUv).xyz;"
        "  vec3 n = normalize(texture(uNormalTex, fsUv).xyz);"
        "  vec3 pos = texture(uPositionTex, fsUv).xyz;"

        "  vec3 l = normalize(vec3(-0.7, 0.3, 0.1));"
        "  vec3 v = normalize(uCameraPos - pos);"
        "  vec3 h = normalize(l + v);"

        "  vec3 color ="
        // diffuse
        "  0.7 * albedo.xyz * max(0.0, dot(n.xyz, l)) +"
        // specular
        "  0.4 * pow(max(0.0, dot(h, n)), 32.0) +"
        // ambient.
        "  0.2 * albedo.xyz;"

        "  outColor = vec4(color, 1.0);"
        "}"
        );



    pointLightShader = LoadNormalShader(
        "#version 330\n"

        "layout(location = 0) in vec3 vsPos;"
        "out vec4 fsPos;"

        "uniform mat4 uVp;"
        "uniform float uLightRadius;"
        "uniform vec3 uLightPosition;"

        "void main()"
        "{"
        "  vec4 pos = uVp * vec4((vsPos * uLightRadius) + uLightPosition, 1.0);"

        "  gl_Position = pos;"
        "  fsPos = pos;"
        "}",


        "#version 330\n"

        "uniform sampler2D uColorTex;"
        "uniform sampler2D uNormalTex;"
        "uniform sampler2D uPositionTex;"


        "out vec4 outColor;"

        "in vec4 fsPos;"

        "uniform float uLightRadius;"
        "uniform vec3 uLightPosition;"
        "uniform vec3 uLightColor;"

        "uniform vec3 uCameraPos;"


        "void main() {"

        // get screen-space position of light sphere
        // (remember to do perspective division.)
        "  vec2 uv = (fsPos.xy / fsPos.w) * 0.5 + 0.5;"

        // now we can sample from the gbuffer for every fragment the light sphere covers.
        "  vec3 albedo = texture(uColorTex, uv).xyz;"
        "  vec3 n = normalize(texture(uNormalTex, uv).xyz);"
        "  vec3 pos = texture(uPositionTex, uv).xyz;"

        "  vec3 lightToPosVector = pos.xyz - uLightPosition;"
        "  float lightDist = length(lightToPosVector);"  // position from light.
        "  vec3 l = -lightToPosVector / (lightDist);"

        // implement fake z-test. If too far from light center, then 0.
        "  float ztest = step(0.0, uLightRadius - lightDist);"

        // light attenuation.
        "  float d = lightDist / uLightRadius;"
        "  float attenuation = 1.0 - d;"
        "  vec3 v = normalize(uCameraPos - pos);"
        "  vec3 h = normalize(l + v);"

        "  vec3 color ="
        // diffuse
        "  uLightColor * albedo.xyz * max(0.0, dot(n.xyz, l)) +"
        // specular
        "  uLightColor * 0.4 * pow(max(0.0, dot(h, n)), 12.0); "

        // finally ztest and attenuation.
        "  color *= ztest * attenuation;"

        "  outColor = vec4(color, 1.0);" // done!
        "}"
        );
    // create the gbuffer. first create fbo:
    GL_C(glGenFramebuffers(1, &fbo));
    GL_C(glBindFramebuffer(GL_FRAMEBUFFER, fbo));

    // RGBA8 color texture-p
    GL_C(glGenTextures(1, &colorTexture));
    GL_C(glBindTexture(GL_TEXTURE_2D, colorTexture));
    GL_C(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, fbWidth, fbHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
    GL_C(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                GL_TEXTURE_2D, colorTexture, 0));
    GL_C(glBindTexture(GL_TEXTURE_2D, 0));

    //  RGBA16F normal texture.
    GL_C(glGenTextures(1, &normalTexture));
    GL_C(glBindTexture(GL_TEXTURE_2D, normalTexture));
    GL_C(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, fbWidth, fbHeight, 0, GL_RGBA, GL_FLOAT, NULL));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
    GL_C(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,
                                GL_TEXTURE_2D, normalTexture, 0));
    GL_C(glBindTexture(GL_TEXTURE_2D, 0));

    //  RGBA16F position texture.
    GL_C(glGenTextures(1, &positionTexture));
    GL_C(glBindTexture(GL_TEXTURE_2D, positionTexture));
    GL_C(glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, fbWidth, fbHeight, 0, GL_RGBA, GL_FLOAT, NULL));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST));
    GL_C(glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST));
    GL_C(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2,
                                GL_TEXTURE_2D, positionTexture, 0));
    GL_C(glBindTexture(GL_TEXTURE_2D, 0));

    // we need a z-buffer for the gbuffer. but we don't need to read from it.
    // so instead create a renderbuffer.
    GL_C(glGenRenderbuffers(1, &depthRenderbuffer));
    GL_C(glBindRenderbuffer(GL_RENDERBUFFER, depthRenderbuffer));
    GL_C(glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32, fbWidth, fbHeight));
    GL_C(glBindRenderbuffer(GL_RENDERBUFFER, 0));
    GL_C(glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRenderbuffer));

    // specify that we can render to all three attachments.
    // this is very important! It won't work otherwise.
    GLenum tgts[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
    GL_C(glDrawBuffers(3, tgts));

    // make sure nothing went wrong:
    GLenum status;
    GL_C(status = glCheckFramebufferStatus(GL_FRAMEBUFFER));
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        printf("Framebuffer not complete. Status: %d", status);
        exit(1);
    }
    GL_C(glBindFramebuffer(GL_FRAMEBUFFER, 0));
    // gbuffer done!

    LoadModel(); // load the obj model.
    CreateSphere(); // create the light sphere geometry.

    while (!glfwWindowShouldClose(window)) {
        float frameStartTime = (float)glfwGetTime();

        glfwPollEvents();
        Render();
        HandleInput();
        glfwSwapBuffers(window);

        //
        // frame rate regulation code:
        //
        float frameEndTime = (float)glfwGetTime();
        float frameDuration = frameEndTime - frameStartTime;
        float sleepDuration = 1.0f / (float)FRAME_RATE - frameDuration;
        if (sleepDuration > 0.0)
            std::this_thread::sleep_for(std::chrono::milliseconds((int)sleepDuration));
        totalTime += 1.0f / (float)FRAME_RATE;
    }

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
