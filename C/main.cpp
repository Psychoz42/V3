#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <optional>
#include <functional>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>


double matrixRes2(double n1_1, double n1_2, double n2_1, double n2_2);
double matrixRes3(double n1_1, double n1_2, double n1_3, double n2_1, double n2_2, double n2_3, double n3_1, double n3_2, double n3_3);
void errorEnd(int err);
bool isEqual(double a, double b);
bool isEqualD(double a, double b);
float radToGrad(float rad);

// Classes start ----------------------------------

class Position{
    public:
        double x;
        double y;
        double z;

        Position(double x_s, double y_s, double z_s){
            x = x_s;
            y = y_s;
            z = z_s;
        }
        Position(){
            x = 0;
            y = 0;
            z = 0;
        }
        std::string text(){
            return std::string{std::to_string(x) + " | " + std::to_string(y) + " | " + std::to_string(z)};
        }
};

class Vector{
    public:
        double vx;
        double vy;
        double vz;

        Vector(double vx_s, double vy_s, double vz_s){
            vx = vx_s;
            vy = vy_s;
            vz = vz_s;
        }
        Vector(){
            vx = 0;
            vy = 0;
            vz = 0;
        }
        std::string text(){
            return std::string{std::to_string(vx) + " | " + std::to_string(vy) + " | " + std::to_string(vz)};
        }
};

double lenV(const Vector &vec);
Vector makeVec(const Position &pos1, const Position &pos2);
double vecsDeg(const Vector &vec1, const Vector &vec2);
double dist(const Position &pos1, const Position &pos2);
Vector vecToUvec(const Vector &vec);
void polsBordering(Position &camPos, Vector &camVec, Vector &pre, int mH, int mW, Vector &xyPerV, Vector &zPerV, float WHK);
Vector precalcVec(Vector &vec, Vector &xyPerV, Vector &zPerV, float WHK);

struct Color
{
    int r, g, b;
};

/*
class polTexture{
    public:
        int polIndex;
        Position anchorDot;
        Vector anchorVec;
        int size;
        std::vector<std::vector<Color>> textureArr;

        polTexture(int polIndex_s, Position anchorDot_s, Vector anchorVec_s, int size_s){
            polIndex = polIndex_s;
            anchorDot = anchorDot_s;
            anchorVec = anchorVec_s;
            size = size_s;
        }

        void addTVec(std::vector<std::vector<Color>> &tArr){
            textureArr = tArr;
        }
};
*/

class Light{
    public:
        Position pos;
        float brightness;
        bool isDir;
        Vector dirVec;

        Light(const Position &pos_s, float brightness_s, bool isDir_s = false, Vector dirVec_s = Vector(0, 0, -1)){
            pos = pos_s;
            brightness = brightness_s;
            isDir = isDir_s;
            dirVec = dirVec_s;
        }
};

class Ray{
    public:
        Position pos;
        Vector vec;

        Ray(const Position &pos_s, const Vector &vec_s){
            pos = pos_s;
            vec = vec_s;
        }
    Ray() = default;
};

void rayToSceneCollusion(const Ray& curRay, Color& col,
    std::string* finalColObj = nullptr,
    Position* finalColPos = nullptr,
    std::vector<std::string>* ignoredObjects = nullptr);

class Vertex{
    public:
        Position pos;

        Vertex(const Position &pos_s){
            pos = pos_s;
        }
        Vertex(){
            pos = Position();
        }
};

class Plane{
    public:
        Vertex v1;
        Vertex v2;
        Vertex v3;

        Vector normalV;

        float free;

        Plane(const Vertex &v1_s, const Vertex &v2_s, const Vertex &v3_s){
            v1 = v1_s;
            v2 = v2_s;
            v3 = v3_s;


            normalV.vx = matrixRes2(v2.pos.y - v1.pos.y, v3.pos.y - v1.pos.y, v2.pos.z - v1.pos.z, v3.pos.z - v1.pos.z);
            normalV.vy = -1 * matrixRes2(v2.pos.x - v1.pos.x, v3.pos.x - v1.pos.x, v2.pos.z - v1.pos.z, v3.pos.z - v1.pos.z);
            normalV.vz = matrixRes2(v2.pos.x - v1.pos.x, v3.pos.x - v1.pos.x, v2.pos.y - v1.pos.y, v3.pos.y - v1.pos.y);

            free = 0;
            free += normalV.vx * -1*v1.pos.x;
            free += normalV.vy * -1*v1.pos.y;
            free += normalV.vz * -1*v1.pos.z;
        }

        Plane(){
            v1 = Vertex();
            v2 = Vertex();
            v3 = Vertex();
        }

        bool haveDot(const Position &pos){
            double pl = matrixRes3(
            pos.x - v1.pos.x, v2.pos.x - v1.pos.x, v3.pos.x - v1.pos.x,
            pos.y - v1.pos.y, v2.pos.y - v1.pos.y, v3.pos.y - v1.pos.y,
            pos.z - v1.pos.z, v2.pos.z - v1.pos.z, v3.pos.z - v1.pos.z);

            if (isEqual(pl, 0.0))
                return true;
            else
                return false;
        }        
};


class Polygon{
    public:
        Color color{0, 0, 0};
        Vertex v1;
        Vertex v2;
        Vertex v3;
        float reflectionK;
        bool haveTexture;
        double minXForVec;
        double maxXForVec;
        double minYForVec;
        double maxYForVec;
        double minZForVec;
        double maxZForVec;

        Plane pl;

        Polygon(const Vertex &v1_s, const Vertex &v2_s, const Vertex &v3_s, const Color &color_s, float rK){ //, bool haveTexture_s
            v1 = v1_s;
            v2 = v2_s;
            v3 = v3_s;

            pl = Plane(v1, v2, v3);
            color = color_s;
            reflectionK = rK;
            //haveTexture = haveTexture_s;
        }
    
        bool haveDot(const Position &pos){
            if (!pl.haveDot(pos))
                return false;
            Vector toDotVec1 = makeVec(v1.pos, pos);
            Vector toDotVec2 = makeVec(v2.pos, pos);
            Vector v1Line1Vec = makeVec(v1.pos, v2.pos);
            Vector v1Line2Vec = makeVec(v1.pos, v3.pos);
            Vector v2Line1Vec = makeVec(v2.pos, v1.pos);
            Vector v2Line2Vec = makeVec(v2.pos, v3.pos);

            if (vecsDeg(toDotVec1, v1Line2Vec) <= vecsDeg(v1Line1Vec, v1Line2Vec))
                if (vecsDeg(toDotVec1, v1Line1Vec) <= vecsDeg(v1Line2Vec, v1Line1Vec))
                    if (vecsDeg(toDotVec2, v2Line2Vec) <= vecsDeg(v2Line1Vec, v2Line2Vec))
                        if (vecsDeg(toDotVec2, v2Line1Vec) <= vecsDeg(v2Line2Vec, v2Line1Vec))
                            return true;
            return false;
        }

        bool haveVec(const Vector &vec){
            if (minXForVec <= vec.vx && vec.vx <= maxXForVec)
                if (minYForVec <= vec.vy && vec.vy <= maxYForVec)
                    if (minZForVec <= vec.vz && vec.vz <= maxZForVec)
                        return true;
            return false;
        }
};

std::optional<Position> rayPolCollusion(const Ray &ray, const Polygon &pol);

std::vector<Polygon> polygons;
std::vector<Light> lights;
//std::vector<polTexture> polTextures;

class Sphere{
    public:
        Position pos;
        float size;
        Color color{0, 0, 0};

        Sphere(const Position &pos_s, float size_s, const Color &color_s){
            pos = pos_s;
            size = size_s;
            color = color_s;
        }
};

std::optional<Position> raySphereCollusion(const Ray &ray, const Sphere &sph);
std::vector<Sphere> spheres;

std::vector<std::string> renderable;

class Camera{
    public:
        Position pos;
        Vector vec;
        int matrixW;
        int matrixH;
        float WHK;
        float FOV;
        float FOVk;
        std::vector<std::string> ignorances{};
        Vector precalcRayVec;
        Vector xyPerV;
        Vector zPerV;
        bool doPostProcessing;
        //Render vars

        Camera(const Position &pos_s, const Vector &vec_s, int mW, int mH, float FOV_s, bool doPostProcessing_s){
            pos = pos_s;
            vec = vecToUvec(vec_s);
            matrixW = mW;
            matrixH = mH;
            WHK = float(mW)/mH;
            FOV = FOV_s;
            FOVk = std::tan(FOV_s * M_PI / 360);
            doPostProcessing = doPostProcessing_s;
        }

        struct Task {
            int x, y;
            Ray curRay;
            Task() = default;
            Task(int x_, int y_, const Ray &ray) : x(x_), y(y_), curRay(ray) {}
        };



        void pixelRender(const Ray &curRay, Color &col){
            std::string finObj;
            Position finPos;
            rayToSceneCollusion(curRay, col, &finObj, &finPos, &ignorances);
            //std::cerr << finObj << std::endl;
            //std::cerr << finPos.text() << std::endl;
        }

        void render(){
            //Prepairing
            std::vector<std::vector<Color>> imageArr(matrixH, std::vector<Color>(matrixW, Color{0,0,0}));

            if(vec.vx != 0 || vec.vy != 0){
                Plane xyPerPlane(Vertex(pos), Vertex(Position(pos.x, pos.y, pos.z + 1.0)), Vertex(Position(pos.x + vec.vx, pos.y + vec.vy, pos.z + vec.vz)));
                xyPerV = vecToUvec(xyPerPlane.normalV);
                Plane zPerPlane(Vertex(pos), Vertex(Position(pos.x + vec.vx, pos.y + vec.vy, pos.z + vec.vz)), Vertex(Position(pos.x + xyPerV.vx, pos.y + xyPerV.vy, pos.z)));
                zPerV = vecToUvec(zPerPlane.normalV);
            }else
                errorEnd(1);

            int numThreads = std::thread::hardware_concurrency();
            std::queue<Task> tasks;
            std::mutex queueMutex;
            std::condition_variable cv;
            bool done = false;

            precalcRayVec = Vector(vec.vx + (-0.5)*xyPerV.vx*WHK, vec.vy + (-0.5)*xyPerV.vy*WHK, vec.vz + 0.5*zPerV.vz);

            //polsBordering(pos, vec, precalcRayVec, matrixH, matrixW, xyPerV, zPerV, WHK);

            //Start of rendering
            for(int y = 0; y < matrixH; y++){
                for(int x = 0; x < matrixW; x++){

                    float curRayVx = precalcRayVec.vx + x/float(matrixW) * xyPerV.vx*WHK + 0.5*zPerV.vx - y/float(matrixH) * zPerV.vx;
                    float curRayVy = precalcRayVec.vy + x/float(matrixW) * xyPerV.vy*WHK + 0.5*zPerV.vy - y/float(matrixH) * zPerV.vy;
                    float curRayVz = precalcRayVec.vz - y/float(matrixH) * zPerV.vz;

                    Ray curRay(pos, vecToUvec(Vector(curRayVx, curRayVy, curRayVz)));
                    
                    tasks.push({ x, y, curRay });
                }
                //std::cerr << round((y*100/float(matrixH)) * 100) / 100.0 << '%' << std::endl;
            }

            std::vector<std::thread> workers;

            auto worker = [&]() {
                while (true) {
                    Task task;
                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        cv.wait(lock, [&]() { return !tasks.empty() || done; });

                        if (done && tasks.empty()) return;

                        task = tasks.front();
                        tasks.pop();
                    }
                    pixelRender(task.curRay, imageArr[task.y][task.x]);
                }
            };

            for (int i = 0; i < numThreads; ++i) {
                workers.emplace_back(worker);
            }

            {
                std::unique_lock<std::mutex> lock(queueMutex);
                done = true;
            }
            cv.notify_all();

            for (auto &t : workers) {
                t.join();
            }
            
            for(int y = 0; y < matrixH; y++){
                for(int x = 0; x < matrixW; x++){
                    std::cout << imageArr[y][x].r << '|' << imageArr[y][x].g << '|' << imageArr[y][x].b << std::endl;
                }
                std::cout << '-' << std::endl;
            }
            
        }
};

// Classes end ------------------------------------

// Funcs start ------------------------------------

bool isEqual(double a, double b){
    //float epsilon = 1e-5;
    return std::fabs(a - b) < 1e-6;
}

double lenV(const Vector &vec){
    return std::sqrt(vec.vx * vec.vx + vec.vy * vec.vy + vec.vz * vec.vz);
}

double matrixRes2(double n1_1, double n1_2, double n2_1, double n2_2){
    return n1_1 * n2_2 - n1_2 * n2_1;
}

double matrixRes3(double n1_1, double n1_2, double n1_3, double n2_1, double n2_2, double n2_3, double n3_1, double n3_2, double n3_3){
    return n1_1 * matrixRes2(n2_2, n2_3, n3_2, n3_3) - n2_1 * matrixRes2(n1_2, n1_3, n3_2, n3_3) + n3_1 * matrixRes2(n1_2, n1_3, n2_2, n2_3);
}

double skal(const Vector &vec1, const Vector &vec2){
    return vec1.vx * vec2.vx + vec1.vy * vec2.vy + vec1.vz * vec2.vz;
}

Vector makeVec(const Position &pos1, const Position &pos2){
    return Vector(pos2.x - pos1.x, pos2.y - pos1.y, pos2.z - pos1.z);
}

double vecsDeg(const Vector &vec1, const Vector &vec2){
    return std::acos(skal(vec1, vec2) / (lenV(vec1) * lenV(vec2)));
}

double dist(const Position &pos1, const Position &pos2){
    return std::sqrt((pos2.x - pos1.x) * (pos2.x - pos1.x) + (pos2.y - pos1.y) * (pos2.y - pos1.y) + (pos2.z - pos1.z) * (pos2.z - pos1.z));
}

Vector vecToUvec(const Vector &vec){
    double vecL = lenV(vec);
    return Vector(vec.vx / vecL, vec.vy / vecL, vec.vz / vecL);
}

Vector precalcVec(Vector &vec, Vector &xyPerV, Vector &zPerV, float WHK){
    return Vector(vec.vx + (-0.5)*xyPerV.vx*WHK, vec.vy + (-0.5)*xyPerV.vy*WHK, vec.vz + 0.5*zPerV.vz);
}



void polsBordering(Position &camPos, Vector &camVec, Vector &pre, int mH, int mW, Vector &xyPerV, Vector &zPerV, float WHK){
    for(std::string r : renderable){

        Polygon *i{nullptr};
        bool toProcess = false;
        
        if (r.substr(0,1) == "p"){
            i = &polygons[std::stoi(r.substr(1, r.size()-1))];
            toProcess = true;
        }
        if (toProcess){
            Vector vecToV1 = vecToUvec(makeVec(camPos, i->v1.pos));
            Vector vecToV2 = vecToUvec(makeVec(camPos, i->v2.pos));
            Vector vecToV3 = vecToUvec(makeVec(camPos, i->v3.pos));

            std::cerr << vecToV1.text() << std::endl;
            std::cerr << vecToV2.text() << std::endl;
            std::cerr << vecToV3.text() << std::endl;
            //std::cerr << camVec.text() << std::endl;
            int y1 = round(((pre.vz - vecToV1.vz) * mH) / zPerV.vz);
            int y2 = round(((pre.vz - vecToV2.vz) * mH) / zPerV.vz);
            int y3 = round(((pre.vz - vecToV3.vz) * mH) / zPerV.vz);
            int x1 = round(((vecToV1.vx - pre.vx - (0.5 - y1/float(mH)) * zPerV.vx) * mW) / (xyPerV.vx * WHK));
            int x2 = round(((vecToV2.vx - pre.vx - (0.5 - y2/float(mH)) * zPerV.vx) * mW) / (xyPerV.vx * WHK));
            int x3 = round(((vecToV3.vx - pre.vx - (0.5 - y3/float(mH)) * zPerV.vx) * mW) / (xyPerV.vx * WHK));
            std::cerr << x1 << " | " << y1 << std::endl;
            std::cerr << x2 << " | " << y2 << std::endl;
            std::cerr << x3 << " | " << y3 << std::endl;
            /*
            i->minXForVec = vecToV1.vx;
            if (vecToV2.vx < i->minXForVec)
                i->minXForVec = vecToV2.vx;
            if (vecToV3.vx < i->minXForVec)
                i->minXForVec = vecToV3.vx;

            i->maxXForVec = vecToV1.vx;
            if (vecToV2.vx > i->maxXForVec)
                i->maxXForVec = vecToV2.vx;
            if (vecToV3.vx > i->maxXForVec)
                i->maxXForVec = vecToV3.vx;
            //
            i->minYForVec = vecToV1.vy;
            if (vecToV2.vy < i->minYForVec)
                i->minYForVec = vecToV2.vy;
            if (vecToV3.vy < i->minYForVec)
                i->minYForVec = vecToV3.vy;

            i->maxYForVec = vecToV1.vy;
            if (vecToV2.vy > i->maxYForVec)
                i->maxYForVec = vecToV2.vy;
            if (vecToV3.vy > i->maxYForVec)
                i->maxYForVec = vecToV3.vy;
            //
            i->minZForVec = vecToV1.vz;
            if (vecToV2.vz < i->minZForVec)
                i->minZForVec = vecToV2.vz;
            if (vecToV3.vz < i->minZForVec)
                i->minZForVec = vecToV3.vz;

            i->maxZForVec = vecToV1.vz;
            if (vecToV2.vz > i->maxZForVec)
                i->maxZForVec = vecToV2.vz;
            if (vecToV3.vz > i->maxZForVec)
                i->maxZForVec = vecToV3.vz;
            

            std::cerr << i->minXForVec << " | " << i->maxXForVec << std::endl;
            std::cerr << i->minYForVec << " | " << i->maxYForVec << std::endl;
            std::cerr << i->minZForVec << " | " << i->maxZForVec << std::endl;
            std::cerr << camVec.text() << std::endl;
            */
        }
    }
}

//float radToGrad(float rad){
//    return (rad / M_PI) * 180;
//}

std::optional<Position> rayPolCollusion(const Ray &ray, const Polygon &pol){
    double t = ray.vec.vx * pol.pl.normalV.vx + ray.vec.vy * pol.pl.normalV.vy + ray.vec.vz * pol.pl.normalV.vz;
    double f = ray.pos.x * pol.pl.normalV.vx + ray.pos.y * pol.pl.normalV.vy + ray.pos.z * pol.pl.normalV.vz + pol.pl.free;

    double t0 = (-1 * f) / t;

    if(t0 < 0.0)
        return std::nullopt;

    double x = ray.vec.vx * t0 + ray.pos.x;
    double y = ray.vec.vy * t0 + ray.pos.y;
    double z = ray.vec.vz * t0 + ray.pos.z;

    if(std::fabs(x - pol.v1.pos.x) > 1e-6 && std::fabs(x - pol.v2.pos.x) > 1e-6 && std::fabs(x - pol.v3.pos.x) > 1e-6)
        if((x > pol.v1.pos.x && x > pol.v2.pos.x && x > pol.v3.pos.x) || (x < pol.v1.pos.x && x < pol.v2.pos.x && x < pol.v3.pos.x))
            return std::nullopt;
    if(std::fabs(y - pol.v1.pos.y) > 1e-6 && std::fabs(y - pol.v2.pos.y) > 1e-6 && std::fabs(y - pol.v3.pos.y) > 1e-6)
        if((y > pol.v1.pos.y && y > pol.v2.pos.y && y > pol.v3.pos.y) || (y < pol.v1.pos.y && y < pol.v2.pos.y && y < pol.v3.pos.y))
            return std::nullopt;
    if(std::fabs(z - pol.v1.pos.z) > 1e-6 && std::fabs(z - pol.v2.pos.z) > 1e-6 && std::fabs(z - pol.v3.pos.z) > 1e-6)
        if((z > pol.v1.pos.z && z > pol.v2.pos.z && z > pol.v3.pos.z) || (z < pol.v1.pos.z && z < pol.v2.pos.z && z < pol.v3.pos.z))
            return std::nullopt;
    
    return Position(x, y, z);
}

std::optional<Position> raySphereCollusion(const Ray &ray, const Sphere &sph){
    double a = ray.vec.vx * ray.vec.vx + ray.vec.vy * ray.vec.vy + ray.vec.vz * ray.vec.vz;
    double b = 2*(ray.pos.x - sph.pos.x)*ray.vec.vx + 2*(ray.pos.y - sph.pos.y)*ray.vec.vy + 2*(ray.pos.z - sph.pos.z)*ray.vec.vz;
    double c = (ray.pos.x - sph.pos.x) * (ray.pos.x - sph.pos.x) + (ray.pos.y - sph.pos.y) * (ray.pos.y - sph.pos.y) + (ray.pos.z - sph.pos.z) * (ray.pos.z - sph.pos.z) - sph.size * sph.size;

    double D = b * b - 4*a*c;
    if (D < 0)
        return std::nullopt;

    if (D == 0){
        double t = (-1*b) / (2*a);
        if (t < 0.0)
            return std::nullopt;
        double x = ray.vec.vx * t + ray.pos.x;
        double y = ray.vec.vy * t + ray.pos.y;
        double z = ray.vec.vz * t + ray.pos.z;
        return Position(x, y, z);
    }

    double t1 = (-1*b - std::sqrt(D)) / (2 * a);
    double t2 = (-1*b + std::sqrt(D)) / (2 * a);

    double x1 = ray.vec.vx * t1 + ray.pos.x;
    double x2 = ray.vec.vx * t2 + ray.pos.x;
    double y1 = ray.vec.vy * t1 + ray.pos.y;
    double y2 = ray.vec.vy * t2 + ray.pos.y;
    double z1 = ray.vec.vz * t1 + ray.pos.z;
    double z2 = ray.vec.vz * t2 + ray.pos.z;

    if(dist(Position(x1, y1, z1), ray.pos) < dist(Position(x2, y2, z2), ray.pos))
        if (t1 > 0.0)
                return Position(x1, y1, z1);
    else
        if (t2 > 0.0)
                return Position(x2, y2, z2);
    
    return std::nullopt;
}

void rayToSceneCollusion(const Ray& curRay, Color& col, std::string* finalColObj, Position* finalColPos, std::vector<std::string>* ignoredObjects)
    {
    float minDist = std::numeric_limits<float>::max();
    int finalRenderObjType;
    Color finalObjColor;
    Position colDot;
    bool colFound = false;
    float light = 0.1;
    std::string fr;

    for(std::string r : renderable){

        Polygon *i{nullptr};
        Sphere *s{nullptr};
        int renderObjType;
        Position tempColDot;
        
        if (r.substr(0,1) == "p"){
            i = &polygons[std::stoi(r.substr(1, r.size()-1))];
            renderObjType = 1;
        }else if (r.substr(0,1) == "s"){
            s = &spheres[std::stoi(r.substr(1, r.size()-1))];
            renderObjType = 2;
        }

        if (renderObjType == 1){
            //if(i->haveVec(curRay.vec)){
                if(auto coordsP = rayPolCollusion(curRay, *i))
                tempColDot = *coordsP;
                else
                    continue;

                if(!i->haveDot(tempColDot))
                    continue;
            //}else
            //    continue;
        }else if (renderObjType == 2){
            if(auto coordsS = raySphereCollusion(curRay, *s))
                tempColDot = *coordsS;
            else
                continue;
        }

        bool skipFlag = false;
        
        if (ignoredObjects){
            for(const std::string ign : *ignoredObjects)
                if(r == ign){
                    skipFlag = true;
                    break;
                }
            
            if(skipFlag)
                continue;
        }

        float newDist = dist(curRay.pos, tempColDot);
        if(newDist <= minDist){
            minDist = newDist;
            colDot = tempColDot;
            fr = r;
            finalRenderObjType = renderObjType;
            if(renderObjType == 1)
                finalObjColor = i->color;
            else
                finalObjColor = s->color;
        }else
            continue;
        
        //if(isEqual(tempColDot.x, curRay.pos.x) && isEqual(tempColDot.y, curRay.pos.y) && isEqual(tempColDot.z, curRay.pos.z)){
        //    continue;
        //}
        colFound = true;
    }
    if(finalColObj)
        *finalColObj = fr;
    
    if(finalColPos)
        *finalColPos = colDot;

    if(colFound){
        for(Light j : lights)
            if (!j.isDir){
                bool haveCol = false;

                for(std::string k : renderable){
                    Polygon *lObjP{nullptr};
                    Sphere *lObjS{nullptr};

                    int lRenderObjType;

                    if (k.substr(0,1) == "p"){
                        lObjP = &polygons[std::stoi(k.substr(1, k.size()-1))];
                        lRenderObjType = 1;
                    }else if (k.substr(0,1) == "s"){
                        lObjS = &spheres[std::stoi(k.substr(1, k.size()-1))];
                        lRenderObjType = 2;
                    }
                    Position lColDot;

                    if (lRenderObjType == 1){
                        if(auto lCoordsP = rayPolCollusion(Ray(j.pos, makeVec(j.pos, colDot)), *lObjP))
                            lColDot = *lCoordsP;
                        else
                            continue;

                        if(!lObjP->haveDot(lColDot))
                            continue;

                    }else if (lRenderObjType == 2){
                        if(auto lCoordsS = raySphereCollusion(Ray(j.pos, makeVec(j.pos, colDot)), *lObjS))
                            lColDot = *lCoordsS;
                        else
                            continue;
                    }

                    if(isEqual(lColDot.x, colDot.x) && isEqual(lColDot.y, colDot.y) && isEqual(lColDot.z, colDot.z))
                        continue;

                    if(dist(j.pos, colDot) > dist(j.pos, lColDot)){
                        haveCol = true;
                        break;
                    }

                }

                if(!haveCol)
                    light += j.brightness / std::pow(dist(j.pos, colDot), 2);
            }

        if(light > 1.0)
            light = 1.0;

        if(finalRenderObjType == 1){
            col.r = round(finalObjColor.r * light);
            col.g = round(finalObjColor.g * light);
            col.b = round(finalObjColor.b * light);
        }else{
            col.r = round(finalObjColor.r * light);
            col.g = round(finalObjColor.g * light);
            col.b = round(finalObjColor.b * light);
        }
    }
}

std::vector<std::string> stringSplit(const std::string &str, const char &diviner){
    std::vector<std::string> splittedStrings;
    int oldSpacePos = 0;
    int spacePos = str.find(diviner);

    while(spacePos != std::string::npos){
        splittedStrings.push_back(str.substr(oldSpacePos, spacePos - oldSpacePos));
        oldSpacePos = spacePos + 1;
        spacePos = str.find(diviner, oldSpacePos);
    }
    splittedStrings.push_back(str.substr(oldSpacePos, spacePos - oldSpacePos));
    return splittedStrings;
}

void errorEnd(int err){
    switch (err)
    {
    case 1:
        std::cerr << "Error: wrong camera vector. X and Y components of vector equals 0. (1)" << std::endl;
        break;
    }

    
    std::exit(1);
}

Camera unloadAll(){

    std::string line;
    int readType = 0;
    bool doPostProcessing_s = false;
    
    while (std::getline(std::cin, line))
    {
        //std::cout << line << std::endl;
        if (line == "pols"){
            readType = 1;
            continue;
        }else if (line == "lights"){
            readType = 2;
            continue;
        }else if (line == "sphs"){
            readType = 3;
            continue;
        //}else if (line == "polT"){
        //    readType = 4;
        //    continue;
        }else if (line == "cam"){
            readType = 4;
            continue;
        }else if (line == "post"){
            readType = 5;
            continue;
        }else{
            int len = line.length();
            std::vector<std::string> lineData = stringSplit(line, ' ');
            switch (readType)
            {
                case 1:
                    {
                    Vertex v1(Position(std::stod(lineData[0]), std::stod(lineData[1]), std::stod(lineData[2])));
                    Vertex v2(Position(std::stod(lineData[3]), std::stod(lineData[4]), std::stod(lineData[5])));
                    Vertex v3(Position(std::stod(lineData[6]), std::stod(lineData[7]), std::stod(lineData[8])));
                    Color color{std::stoi(lineData[9]), std::stoi(lineData[10]), std::stoi(lineData[11])};
                    polygons.push_back(Polygon(v1, v2, v3, color, std::stod(lineData[12]))); //, std::stoi(lineData[13])
                    }
                    break;
                case 2:
                    {
                    Position pos(std::stod(lineData[0]), std::stod(lineData[1]), std::stod(lineData[2]));
                    Vector dirVec(std::stod(lineData[5]), std::stod(lineData[6]), std::stod(lineData[7]));
                    lights.push_back(Light(pos, std::stod(lineData[3]), std::stoi(lineData[4]), dirVec));
                    }
                    break;
                case 3:
                    {
                    Position pos(std::stod(lineData[0]), std::stod(lineData[1]), std::stod(lineData[2]));
                    Color color{std::stoi(lineData[4]), std::stoi(lineData[5]), std::stoi(lineData[6])};
                    spheres.push_back(Sphere(pos, std::stod(lineData[3]), color));
                    }
                    break;
                /*
                case 4:
                    {
                    Position pos(std::stod(lineData[1]), std::stod(lineData[2]), std::stod(lineData[3]));
                    Vector vec(std::stod(lineData[4]), std::stod(lineData[5]), std::stod(lineData[6]));
                    polTextures.push_back(polTexture(std::stoi(lineData[0]), pos, vec, std::stoi(lineData[7])));
                    readType = 6;
                    }
                    break;
                */
                case 4: //5
                    {
                    Position pos(std::stod(lineData[0]), std::stod(lineData[1]), std::stod(lineData[2]));
                    Vector vec(std::stod(lineData[3]), std::stod(lineData[4]), std::stod(lineData[5]));

                    for(int i = 0; i < polygons.size(); i++){
                        std::string s{"p"};
                        s.append(std::to_string(i));
                        renderable.push_back(s);
                    }
                    for (int i = 0; i < spheres.size(); i++){
                        std::string s{"s"};
                        s.append(std::to_string(i));
                        renderable.push_back(s);
                    }
                    return Camera(pos, vec, std::stoi(lineData[6]), std::stoi(lineData[7]), std::stod(lineData[8]), doPostProcessing_s);
                    }
                    break;
                case 5:
                    {
                        doPostProcessing_s = true;
                    }
                    break;
                /*
                case 6:
                    {
                    std::vector<std::vector<Color>> tArr;
                    std::vector<Color> tempArr;
                    Color tempCol;
                    int i;
                    for(std::string str : lineData){
                        if(str == "|"){
                            tArr.push_back(tempArr);
                            tempArr.clear();
                            i = 0;
                        }
                        if(i < 3){
                            switch (i)
                            {
                            case 0:
                                {
                                tempCol.r = std::stoi(str);
                                }
                                break;
                            case 1:
                                {
                                tempCol.g = std::stoi(str);
                                }
                                break;
                            case 2:
                                {
                                tempCol.b = std::stoi(str);
                                }
                                break;
                            }
                        } else {
                            tempArr.push_back(tempCol);
                        }
                    }
                    polTextures[-1].addTVec(tArr);
                    }
                    break;
                */
            }
        }
    }
}

// Funcs end --------------------------------------

int main() {
    
    Camera cam = unloadAll();

    cam.render();

}