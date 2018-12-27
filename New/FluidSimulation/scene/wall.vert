# version 440

layout (location = 0) in vec3 v;

uniform vec3 Eye;

uniform mat4 MVP;

out vec3 Color;

void main()
{
    Color = vec3(0.8, 0.4, 0.0);
    gl_Position = MVP * vec4(v, 1.0);
}
