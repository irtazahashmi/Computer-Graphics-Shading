#pragma once
// Disable warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()
#include <iostream>

// !!! DO NOT MODIFY THIS STRUCT !!!
struct ShadingData {
    glm::vec3 Kd { 0.5f, 0.5f, 0.5f }; // Diffuse coefficient per vertex.
    glm::vec3 Ks { 0.5f, 0.5f, 0.5f }; // Specularity coefficient per vertex.
    float shininess = 20.0f; // Exponent for Phong and Blinn-Phong specularities per vertex.
    int toonDiscretize = 4; // Number of levels in Toon shading.
    float toonSpecularThreshold = 0.49f; // Threshold for specularity.
};

// ==========================
// =====    EXERCISE    =====
// ==========================

// For debugging purposes or variable changes (e.g., modify the Toon threshold as done below).
// Please notice that some keys are already in use!
void yourKeyboardFunction(unsigned char key, 
                          ShadingData& shadingData)
{
    std::cout << "Key not used so far, so executing your code!" << std::endl;

    //recolor the mesh
    switch (key) {
    case 't':
        shadingData.toonSpecularThreshold -= 0.001f;
        break;
    case 'T':
        shadingData.toonSpecularThreshold += 0.001f;
        break;
    case 'd':
        shadingData.toonDiscretize = std::max(1, shadingData.toonDiscretize - 1);
        break;
    case 'D':
        shadingData.toonDiscretize += 1;
        break;
    }

    std::cout << "ToonSpecular" << shadingData.toonSpecularThreshold << std::endl;
    std::cout << "Toon Discretization levels" << shadingData.toonDiscretize << std::endl;
}

// Debug function.
glm::vec3 debugColor(const ShadingData& shadingData, 
                        const glm::vec3& vertexPos, 
                        const glm::vec3& normal, 
                        const glm::vec3& lightPos, 
                        const glm::vec3& cameraPos)
{
    // This function you can use in any way you like!
    // E.g., for debugging purposes!
    return shadingData.Kd;
}

// Standard lambertian shading: Kd * dot(N,L), clamped to zero when negative. Where L is the light vector.
glm::vec3 diffuseOnly(const ShadingData& shadingData, 
                        const glm::vec3& vertexPos, 
                        const glm::vec3& normal, 
                        const glm::vec3& lightPos)
{

    // Kd - diffuse coefficient
    glm::vec3 kd = shadingData.Kd;

    // L - lightVector = lightPos - vertexPos
    glm::vec3 lightVector = lightPos - vertexPos;
    glm::vec3 lightVectorNormalized = glm::normalize(lightVector);

    // N - normlaized normal
    glm::vec3 normalizedNormal = glm::normalize(normal);

    // dot(N, L)
    float dotProductNormalLightVector = glm::dot(normalizedNormal, lightVectorNormalized);

    
    // If the light source is behind the normal meaning that cos(theta) < 0 -> should result in black 
    if (dotProductNormalLightVector > 0) {
        // using the lambertian formula
        glm::vec3 lambertianModel = kd * dotProductNormalLightVector;
        return lambertianModel;
    } else {
        //if negative, clamped to zero
        return glm::vec3{ 0, 0, 0 };
    }
}

// Phong (!) Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
// Follow the course, only calculate Ks pow(dot(V,R),shininess), where V is the view vector and R is the Reflection 
// vector of the light (like in pool billard computed from the LightPos, vertexPos and normal).
// When computing specularities like this, verify that the light is on the right side of the surface, with respect 
// to the normal
// E.g., for a plane, the light source below the plane cannot cast light on the top, hence, there can also not be any 
// specularity.
glm::vec3 phongSpecularOnly(const ShadingData& shadingData, 
                            const glm::vec3& vertexPos, 
                            const glm::vec3& normal, 
                            const glm::vec3& lightPos, 
                            const glm::vec3& cameraPos)
{
    // specularity coefficient
    glm::vec3 ks = shadingData.Ks;

    // shininess
    float shininess = shadingData.shininess;

    // N - normalized vector
    glm::vec normalizedNormal = glm::normalize(normal);

    // L - finding lightVector = lightPos - vertexPos
    glm::vec3 lightVector = lightPos - vertexPos;
    glm::vec3 lightVectorNormalized = glm::normalize(lightVector);

    // R - Reflection vector of the light = 2 * (L * N) * N - L
    glm::vec3 reflectionVectorOfLight = 2 
        * (glm::dot(lightVectorNormalized, normalizedNormal))
        * normalizedNormal - lightVectorNormalized;

    // V - view vector 
    glm::vec3 viewVector = cameraPos - vertexPos;
    glm::vec3 viewVectorNormalized = glm::normalize(viewVector);

    if (glm::dot(normalizedNormal, lightVectorNormalized) > 0) {
        //Ks pow(dot(V, R), shininess)
        return glm::vec3{ ks * powf(glm::dot(viewVectorNormalized, reflectionVectorOfLight), shininess) };
    }
    else {
        // light placed behind the model should not illuminate it from the front -> black
        return glm::vec3{ 0, 0, 0 };
    }
}

// Blinn-Phong Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
// Be careful!!! The pseudo code does some additional modifications to the formula seen in the course.
// Follow the course version and calculate ONLY Ks * pow(dot(N,H), shininess). The definition of H is 
// given on the page above and in the course.
// The same test as before should be used.
glm::vec3 blinnPhongSpecularOnly(const ShadingData& shadingData, 
                                const glm::vec3& vertexPos, 
                                const glm::vec3& normal, 
                                const glm::vec3& lightPos, 
                                const glm::vec3& cameraPos)
{

    // ks - specularity coefficient
    glm::vec3 ks = shadingData.Ks;

    // shininess
    float shininess = shadingData.shininess;

    // N - normalised normal
    glm::vec3 normalizedNormal = glm::normalize(normal);

    // L - light vector : lightPos - vertexPos
    glm::vec3 lightVector = lightPos - vertexPos;
    glm::vec3 lightVectorNormalized = glm::normalize(lightVector);

    // V - view vector : cameraPos - vertexPos
    glm::vec3 viewVector = cameraPos - vertexPos;
    glm::vec3 viewVectorNormalized = glm::normalize(viewVector);

    // H - half vector: the unit vector exactly between the view direction and the direction towards the light.
    glm::vec3 halfVector = viewVectorNormalized + lightVectorNormalized;
    glm::vec3 halfVectorNormalized = glm::normalize(halfVector);


    if (glm::dot(halfVectorNormalized, normalizedNormal) > 0) {
        //Ks * pow(dot(N,H), shininess)
        return glm::vec3{ ks * powf(glm::dot(halfVectorNormalized, normalizedNormal), shininess) };
    }
    else {
        // light placed behind the model should not illuminate it from the front -> black
        return glm::vec3{0, 0, 0};
    }
}

// Diffuse Toon Shading.
// Use the variable ToonDiscretize.
// Normal diffuse shading has values between 0 and Kd (the diffuse coefficient).
// In toon shading, these values are discretized.
// This means, there are N (ToonDiscretize) uniform intervals between 0 and Kd - in this example, 
// we assume a single color channel, you should use the values from Kd.
// Let c(0)=0, c(1)...c(N) =Kd be the N+1 boundary values of these intervals.
// For a value v in [c(i), c(i+1)), the function should return (c(i)+c(i+1))/2.
// For v=0, return black (0)
// For v=Kd, return Kd
glm::vec3 toonShadingNoSpecular(const ShadingData& shadingData, 
                                const glm::vec3& vertexPos, 
                                const glm::vec3& normal, 
                                const glm::vec3& lightPos)
{
    // kd
    glm::vec3 kd = shadingData.Kd;

    // N - normalized normal
    glm::vec3 normalizedNormal = glm::normalize(normal);

    // L - light vector : lightPos - vertexPos
    glm::vec3 lightVector = lightPos - vertexPos;
    glm::vec3 lightVectorNormalized = glm::normalize(lightVector);

    // dot(N, L)
    float dotProductNormalLightVector = glm::dot(normalizedNormal, lightVectorNormalized);

    // edge case - when v = 0, return 0
    if (dotProductNormalLightVector <= 0) {
        return glm::vec3{ 0,0,0 };
    }

    //For v = Kd, return Kd
    if (dotProductNormalLightVector == kd.x) {
        return kd;
    }

    // n uniform intervals between 0 and kd
    float toonDiscreteValue = (float) shadingData.toonDiscretize;

    // find the right interval according to the threshold
    float threshold = 1 / toonDiscreteValue;

    // get the right bucket
    int i = 0;
    while (dotProductNormalLightVector - threshold > 0) {
        dotProductNormalLightVector -= threshold;
        i++;
    }

    //return (c(i)+c(i+1))/2
    glm::vec3 ci = (i * threshold * kd);
    glm::vec3 ciPlus1 = ((i + 1) * threshold * kd);
    return (ci + ciPlus1) / 2.0f;
}

// Specular Toon shading.
// The goal is to add white highlights. If the Blinn-Phong Specularity (without multiplying by Ks!) 
// has a value larger or equal to ToonSpecularThreshold,
//  then return white (vec3(1)), else return black.
glm::vec3 toonShadingOnlySpecular(const ShadingData& shadingData, 
                                    const glm::vec3& vertexPos, 
                                    const glm::vec3& normal, 
                                    const glm::vec3& lightPos, 
                                    const glm::vec3& cameraPos)
{

    // ks - specularity coefficient
    glm::vec3 ks = shadingData.Ks;

    // shininess
    float shininess = shadingData.shininess;

    // N - normalised normal
    glm::vec3 normalizedNormal = glm::normalize(normal);

    // L - light vector : lightPos - vertexPos
    glm::vec3 lightVector = lightPos - vertexPos;
    glm::vec3 lightVectorNormalized = glm::normalize(lightVector);

    // V - view vector : cameraPos - vertexPos
    glm::vec3 viewVector = cameraPos - vertexPos;
    glm::vec3 viewVectorNormalized = glm::normalize(viewVector);

    // H - half vector: the unit vector exactly between the view direction and the direction towards the light.
    glm::vec3 halfVector = viewVectorNormalized + lightVectorNormalized;
    glm::vec3 halfVectorNormalized = glm::normalize(halfVector);
   
    float blinnPhongSpecularity = powf(glm::dot(halfVectorNormalized, normalizedNormal), shininess);

    //  If blinnPhongSpecularity has a value larger or equal to ToonSpecularThreshold, then return white (vec3(1))
    if (blinnPhongSpecularity >= shadingData.toonSpecularThreshold) {
        return glm::vec3{ 1, 1, 1 };
    } else {
        //else return black.
        return glm::vec3{ 0, 0, 0 };
    }
}

// RETURN the new light position, defined as follows:
// selectedPos is a location on the mesh. Use this location to place the light 
// source to cover the location as seen from camPos (cover the cursor).
// Further, the light should be at a distance of 1.5 from the origin of the scene - 
// in other words, located on a sphere of radius 1.5 around the origin.
// The selectedPos is guaranteed to always be within the sphere.


// helper method for solving quadratic equatuions
float quadraticEquationSolver(const float& a, 
                              const float& b, 
                              const float& c) 
{
    float numerator = -b + sqrtf(powf(b, 2) - 4 * a * c);
    float denominaotor = 2.0f * a;
    return numerator / denominaotor;
}

glm::vec3 userInteractionSphere(const glm::vec3& selectedPos, 
                                const glm::vec3& camPos)
{
    // this article helped me to do this question
    // https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

    // a 
    float a = 1.0f;

    // direction vec = camPos - selectedPos
    glm::vec3 directionVector = camPos - selectedPos;
    glm::vec3 directionVectorNormalized = glm::normalize(directionVector);

    // 2*dot(direction vector, selected pos) -> b
    float dotProductDirectionVectorSelectedPosition = 2 * glm::dot(directionVectorNormalized, selectedPos);
    float b = dotProductDirectionVectorSelectedPosition;

    // distance(camPos-selectedPos)^2 - r^2 -> c
    float radius = 1.5f;
    float distanceToSelectedPos = glm::length(selectedPos);
    float c = powf(distanceToSelectedPos, 2) - powf(radius, 2);

    // finding t 
    float t = quadraticEquationSolver(a, b, c);

    // using the ray equation to find the intersection point: p = o + dt
    glm::vec3 intersectionPoint = selectedPos + directionVector * t;
    return intersectionPoint;
}

// RETURN the new light position such that the light towards the selectedPos is orthgonal to the normal at that location
// --- in this way, the shading boundary will be exactly at this location.
// there are several ways to do this, choose one you deem appropriate given the current light position
// no panic, I will not judge what solution you chose, as long as the above condition is met.
glm::vec3 userInteractionShadow(const glm::vec3& selectedPos, const glm::vec3& selectedNormal, const glm::vec3& lightPos)
{
    return glm::vec3(1, 0, 1);
}

// RETURN the new light position such that a specularity (highlight) will be located at selectedPos, 
// when viewed from cameraPos.
// Please ensure also that the light is at a distance of 1 from selectedPos!
// There is only ONE way of doing this!
glm::vec3 userInteractionSpecular(const glm::vec3& selectedPos, 
                                  const glm::vec3& selectedNormal, 
                                  const glm::vec3& cameraPos)

    // specularity's position depends on the view AND the light
{
    // direction vector
    glm::vec3 directionVector = selectedPos - cameraPos;

    // reflection vector of direction vector
    glm::vec3 reflectionVectorDirection = glm::reflect(directionVector, selectedNormal);
    glm::vec3 reflectionVectorDirectionNormalized = glm::normalize(reflectionVectorDirection);

    // the new ligth position 
    glm::vec3 newLightPosition = reflectionVectorDirectionNormalized + selectedPos;

    // returning the resulting new light position
    return newLightPosition;
}
