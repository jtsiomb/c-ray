uniform mat4 cam_mat;
uniform vec3 cam_pos;
uniform float fov, aspect;
uniform sampler2D geom_tex;
uniform int obj_num, geom_tex_height;

#define RAY_MAG		1000.0

#define TC_PIX0		0.05
#define TC_PIX1		0.3
#define TC_PIX2		0.55

struct ray {
	vec3 orig, dir;
};

struct material {
	vec3 col;
	float spow;
	float refl;
};

struct sphere {
	vec3 pos;
	float rad;
	struct material mat;
};

struct spoint {
	bool hit;
	vec3 pos, norm, vref;
	float dist;
};

vec3 get_sample_pos(float x, float y)
{
	return vec3(aspect * (x * 2.0 - 1.0), y * 2.0 - 1.0, 0.0);
}

ray get_primary_ray(float x, float y)
{
	ray r;
	vec3 dir;
	
	r.orig = vec3(0.0, 0.0, 0.0);
	r.dir = get_sample_pos(x, y);
	r.dir.z = 1.0 / tan(0.5 * fov);
	r.dir *= RAY_MAG;

	r.dir = (cam_mat * vec4(r.dir, 1.0)).xyz;
	r.orig = (cam_mat * vec4(r.orig, 1.0)).xyz + cam_pos;

	r.dir += r.orig;
	return r;
}

#define SQ(x)	((x) * (x))
#define ERR_MARGIN		1e-6

spoint ray_sphere(sphere sph, ray r)
{
	float a, b, c, d, sqrt_d, t1, t2;
	spoint res;

	res.hit = false;

	a = SQ(r.dir.x) + SQ(r.dir.y) + SQ(r.dir.z);
	b = 2.0 * r.dir.x * (r.orig.x - sph.pos.x) +
		2.0 * r.dir.y * (r.orig.y - sph.pos.y) +
		2.0 * r.dir.z * (r.orig.z - sph.pos.z);
	c = SQ(sph.pos.x) + SQ(sph.pos.y) + SQ(sph.pos.z) +
		SQ(r.orig.x) + SQ(r.orig.y) + SQ(r.orig.z) +
		2.0 * (-sph.pos.x * r.orig.x - sph.pos.y * r.orig.y - sph.pos.z * r.orig.z) - SQ(sph.rad);
	d = SQ(b) - 4.0 * a * c;
	if(d < 0.0) {
		return res;
	}

	sqrt_d = sqrt(d);
	t1 = (-b + sqrt_d) / (2.0 * a);
	t2 = (-b - sqrt_d) / (2.0 * a);

	if((t1 < ERR_MARGIN && t2 < ERR_MARGIN) || (t1 > 1.0 && t2 > 1.0)) {
		return res;
	}

	if(t1 < ERR_MARGIN) t1 = t2;
	if(t2 < ERR_MARGIN) t2 = t1;

	res.dist = t1 < t2 ? t1 : t2;
	res.pos = r.orig + r.dir * res.dist;
	res.norm = (res.pos - sph.pos) / sph.rad;
	res.vref = normalize(reflect(r.dir, res.norm));
	res.hit = true;

	return res;
}

void main()
{
	vec4 col = vec4(0.0, 0.0, 0.0, 1.0);

	/*
	sphere sph;
	sph.pos = vec3(0.0, 0.0, 0.0);
	sph.rad = 1.0;
	sph.mat.col = vec3(1.0, 0.0, 0.0);
	sph.mat.spow = 60.0;
	sph.mat.refl = 0.0;
	*/

	ray r = get_primary_ray(gl_TexCoord[0].s, gl_TexCoord[0].t);

	sphere nearest_obj;
	spoint nearest_sp;
	nearest_sp.dist = RAY_MAG;

	for(int i=0; i<obj_num; i++) {
		vec4 pix[4];
		float obj_tc = float(i) / float(geom_tex_height);

		pix[0] = texture2D(geom_tex, vec2(TC_PIX0, obj_tc));
		pix[1] = texture2D(geom_tex, vec2(TC_PIX1, obj_tc));
		pix[2] = texture2D(geom_tex, vec2(TC_PIX2, obj_tc));

		sphere sph;
		sph.pos = pix[0].xyz;
		sph.rad = pix[0].w;
		sph.mat.col = pix[1].xyz;
		sph.mat.spow = pix[2].x;
		sph.mat.refl = pix[2].y;

		spoint sp = ray_sphere(sph, r);
		if(sp.hit && sp.dist < nearest_sp.dist) {
			nearest_sp = sp;
			nearest_obj = sph;
		}
	}

	if(nearest_sp.dist < RAY_MAG) {
		col.xyz = nearest_obj.mat.col;
	}

	gl_FragColor = col;
}
