#version 440

layout (location = 0) in vec3 v;
layout (location = 1) in vec2 vt;
layout (location = 2) in vec3 vn;

out vec3 LightIntensity;

uniform vec3 LightPosition;
uniform vec3 Kd;
uniform vec3 Ld;

uniform mat4 MVP;

void main()
{
    vec3 tnorm = normalize(vn);
    vec3 s = normalize(LightPosition - v);
    LightIntensity = Ld * Kd * max(dot(s, tnorm), 0.0);
    gl_Position = MVP * vec4(v, 1.0);
}
