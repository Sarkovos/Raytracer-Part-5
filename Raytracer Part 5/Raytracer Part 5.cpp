/* Created by: Miguel R. Kunkle
* Date: 2/25/23
* Final step of a simple raytracer.
* Reflections, coming right back at ya.
* Also includes a camera rotation.
*/

#include <SDL.h>
#include <iostream>
#include <cmath>
#include <array>
#include <tuple>
#undef main

//-------------------------------
// Linear Algebra functions
// ------------------------------

//dot product of two 3d vectors
double DotProduct(std::array<double, 3> v1, std::array<double, 3> v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

//computes v1 - v2
std::array<double, 3> subtract(std::array<double, 3> v1, std::array<double, 3> v2)
{
	std::array<double, 3> a;
	for (int i = 0; i < 3; i++)
	{
		a[i] = v1[i] - v2[i];
	}

	return a;
}

//computes v1 + v2
std::array<double, 3> add(std::array<double, 3> v1, std::array<double, 3> v2)
{
	std::array<double, 3> b;
	for (int i = 0; i < 3; i++)
	{
		b[i] = v1[i] + v2[i];
	}

	return b;
}

//length of a 3D vector
double length(std::array<double, 3> vec)
{
	return std::sqrt(DotProduct(vec, vec));
}

//scales a vector by k
std::array<double, 3> multiply(double k, std::array<double, 3> vec)
{
	std::array<double, 3> product;

	for (int i = 0; i < 3; i++)
	{
		product[i] = vec[i] * k;
	}

	return product;
}

//multiplies a matrix and a vector
std::array <double, 3> multiplyMatrix(std::array<std::array<double, 3>, 3> mat, std::array<double, 3> vec)
{
	std::array<double, 3> result = { 0,0,0 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result[i] += vec[j] * mat[i][j];
		}
	}
	return result;
}

// Clamps a color to the canonical color range.
//min makes sure value doesnt go over 255(highest color range), max makes sure value doesnt go below 0
std::array<double, 3> clamp(std::array<double, 3> vec)
{
	std::array<double, 3> clampColor;
	for (int i = 0; i < 3; i++)
	{
		clampColor[i] = std::fmin(255, std::fmax(0, vec[i]));
	}

	return clampColor;
}

// Computes the reflection of v1 respect to v2.
std::array<double, 3> ReflectRay(std::array<double, 3> v1, std::array<double, 3> v2)
{
	return subtract(multiply(2 * DotProduct(v1, v2), v2), v1);
}


//------------------------------------------------------------
// Shiny Raytracer
// -----------------------------------------------------------

//Scene Setup
double viewport_size = 1;
double  projection_plane_z = 1;
std::array<double, 3> camera_position = { 3, 0, 1 };
std::array<double, 3> background_color = { 0, 0, 0 };
double recursion_depth = 3;
double epsilon = 0.001;

std::array<std::array<double, 3>, 3> camera_rotation = { { {0.7071,0,-0.7071}, {0, 1, 0}, {0.7071,0, 0.7071}} };



//canvas width and height
const double WIDTH = 600;
const double HEIGHT = 600;


//--------
//Spheres
//--------

//To create spheres
struct Sphere
{
	std::array<double, 3> center;
	double radius;
	std::array<double, 3> color;
	double specular;
	double reflective;
};

//create spheres
Sphere sphere1{ {0, -1, 3}, 1, {255, 0, 0}, 500, 0.2 };		//Red, shiny, bit reflective
Sphere sphere2{ {2, 0, 4}, 1,{0, 0, 255 }, 500, 0.3 };		//Blue, shiny, more reflective
Sphere sphere3{ {-2, 0, 4}, 1,{0, 255, 0 }, 10, 0.4 };       //Green, slightly shiny, more reflective
Sphere sphere4{ {0,-5001,0}, 5000, {255,0,255}, 1000, 0.5 }; //Purple, very shiny, half reflective

Sphere* spheres[] = { &sphere1, &sphere2, &sphere3, &sphere4 };

//-------
//Lights
//-------

//To create lights
struct Light
{
	std::string ltype;
	double intensity;
	std::array<double, 3> position;
};


//Light array
Light light1{ "AMBIENT", 0.2 };
Light light2{ "POINT", 0.6, {2, 1, 0} };
Light light3{ "DIRECTIONAL", 0.2, {1, 4, 4} };

Light* lights[] = { &light1, &light2, &light3 };




//define "infinity"
const double INFIN = 2147483647;

//Conerts 2D canvas coordinates to 3D viewport coordinates
std::array<double, 3> CanvasToViewport(std::array<double, 2>  p2d)
{
	std::array<double, 3> a;

	a[0] = p2d[0] * viewport_size / WIDTH;
	a[1] = p2d[1] * viewport_size / HEIGHT;
	a[2] = projection_plane_z;

	return a;
}

// Computes the intersection of a ray and a sphere. Returns the values of t for the intersections.
//Camera position is origin
//Direction is what CanvasToViewport returns

std::array<double, 2> IntersectRaySphere(std::array<double, 3> origin, std::array<double, 3> direction, Sphere* sphere)
{
	std::array<double, 3> oc = subtract(origin, sphere->center);

	double k1 = DotProduct(direction, direction);
	double k2 = 2 * DotProduct(oc, direction);
	double k3 = DotProduct(oc, oc) - sphere->radius * sphere->radius;

	double discriminant = k2 * k2 - 4 * k1 * k3;
	if (discriminant < 0)
	{
		std::array<double, 2> infinity = { INFIN, INFIN };
		return infinity;
	}

	double t1 = (-k2 + std::sqrt(discriminant)) / (2 * k1);
	double t2 = (-k2 - std::sqrt(discriminant)) / (2 * k1);

	std::array<double, 2> t = { t1,t2 };

	return t;
}

std::tuple<Sphere*, double> closestIntersection(std::array<double, 3>  origin, std::array<double, 3> direction, double min_t, double max_t)
{
	double closest_t = INFIN;
	Sphere* closest_sphere = nullptr;

	//4 IS HARD CODED FIX LATER LOVE YOU MWAH
	for (int i = 0; i < 4; i++)
	{
		std::array<double, 2> ts = IntersectRaySphere(origin, direction, spheres[i]);

		if (ts[0] < closest_t && min_t < ts[0] && ts[0] < max_t)
		{
			closest_t = ts[0];
			closest_sphere = spheres[i];
		}

		if (ts[1] < closest_t && min_t < ts[1] && ts[1] < max_t)
		{
			closest_t = ts[1];
			closest_sphere = spheres[i];
		}
	}
	if (closest_sphere)
	{
		std::tuple<Sphere*, double> intersection(closest_sphere, closest_t);
		return intersection;
	}

	std::tuple<Sphere*, double> intersectionNull(nullptr, 0);
	return intersectionNull;
}

//Its in the name
double ComputeLighting(std::array<double, 3> point, std::array<double, 3> normal, std::array<double, 3> V, double specular)
{
	double intensity = 0;
	double length_n = length(normal); //Verifying the length of our normal is 1.

	//amount of lights is hard coded to 3
	for (int i = 0; i < 3; i++)
	{
		Light* light = lights[i];


		if (light->ltype == "AMBIENT")
		{
			intensity += light->intensity;
		}

		else
		{
			std::array<double, 3>  vec_l;
			double t_max;
			if (light->ltype == "POINT")
			{
				vec_l = subtract(light->position, point);
				t_max = 1;
			}

			else //This is directional light
			{
				vec_l = light->position;
				t_max = INFIN;
			}

			//Shadow check
			std::tuple<Sphere*, double> blocker = closestIntersection(point, vec_l, epsilon, t_max);
			if (std::get<0>(blocker) != nullptr)
			{
				continue;
			}

			//Diffuse
			double n_dot_l = DotProduct(normal, vec_l);
			if (n_dot_l > 0)
			{
				intensity += light->intensity * n_dot_l / (length_n * length(vec_l));
			}

			//Specular
			if (specular != -1)
			{
				std::array<double, 3> R = ReflectRay(vec_l, normal);
				double r_dot_v = DotProduct(R, V);

				if (r_dot_v > 0)
				{
					intensity += light->intensity * pow(r_dot_v / (length(R) * length(V)), specular);
				}

			}
		}


	}
	return intensity;
}

//Traces a ray against the set of spheres in the scene.
//min_t is 1, max_t is infinity.
std::array<double, 3> TraceRay(std::array<double, 3>  origin, std::array<double, 3> direction, double min_t, double max_t, double depth)
{
	std::tuple<Sphere*, double> intersection = closestIntersection(origin, direction, min_t, max_t);

	if (std::get<0>(intersection) == nullptr)
	{
		return background_color;
	}

	Sphere* closest_sphere = std::get<0>(intersection);
	double closest_t = std::get<1>(intersection);

	std::array<double, 3> point = add(origin, multiply(closest_t, direction));
	std::array<double, 3> normal = subtract(point, closest_sphere->center);
	normal = multiply(1.0 / length(normal), normal);

	std::array<double, 3> view = multiply(-1, direction);
	double lighting = ComputeLighting(point, normal, view, closest_sphere->specular);

	std::array<double, 3>  local_color = multiply(lighting, closest_sphere->color);

	//If we hit the recursion limit or the object is not reflective, we're done
	double r = closest_sphere->reflective;
	if (r <= 0 || depth <= 0)
	{
		return local_color;
	}

	//Compute the reflected color
	std::array<double, 3> reflected_ray = ReflectRay(multiply(-1, direction), normal);
	std::array<double, 3> reflected_color = TraceRay(point, reflected_ray, epsilon, INFIN, depth - 1);



	return add(multiply(1 - r, local_color), multiply(r, reflected_color));
}


void PutPixel(int x, int y, std::array<double, 3> color, SDL_Renderer* renderer)
{
	x = WIDTH / 2 + x;
	y = HEIGHT / 2 - y - 1;

	if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT)
	{
		return;
	}

	SDL_SetRenderDrawColor(renderer, color[0], color[1], color[2], 255);
	SDL_RenderDrawPoint(renderer, x, y);
}

int main()
{
	//SDL Making Window
	SDL_Window* window = nullptr;
	SDL_Renderer* renderer = nullptr;

	SDL_Init(SDL_INIT_VIDEO);
	SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, 0, &window, &renderer);

	SDL_SetRenderDrawColor(renderer, 255, 0, 255, 255);
	SDL_RenderClear(renderer);

	for (double x = -WIDTH / 2; x < WIDTH / 2; x++)
	{
		for (double y = -HEIGHT / 2; y < HEIGHT / 2; y++)
		{
			std::array<double, 2> XY = { x,y };
			std::array<double, 3> direction = CanvasToViewport(XY);
			direction = multiplyMatrix(camera_rotation, direction);
			std::array<double, 3> color = TraceRay(camera_position, direction, 1, INFIN, recursion_depth);
			PutPixel(x, y, clamp(color), renderer);
		}
	}

	SDL_RenderPresent(renderer);

	//delay is currently 3 seconds
	SDL_Delay(3000);
	return 0;
}


