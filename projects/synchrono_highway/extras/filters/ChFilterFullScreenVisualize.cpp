// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Eric Brandt, Asher Elmquist
// =============================================================================
//
// =============================================================================

#include "ChFilterFullScreenVisualize.h"
#include "chrono_sensor/ChOptixSensor.h"
#include "chrono_sensor/utils/CudaMallocHelper.h"

#include <cuda_runtime_api.h>

namespace chrono {
namespace sensor {

int ChFilterFullScreenVisualize::s_windowCount = 0;
std::mutex ChFilterFullScreenVisualize::s_glfwMutex;

CH_SENSOR_API ChFilterFullScreenVisualize::ChFilterFullScreenVisualize(int w, int h, std::string name, bool is_fullscreen) : m_w(w), m_h(h), m_is_fullscreen(is_fullscreen), ChFilter(name) {}

CH_SENSOR_API ChFilterFullScreenVisualize::~ChFilterFullScreenVisualize() {
    ChFilterFullScreenVisualize::OnCloseWindow();
}

CH_SENSOR_API void ChFilterFullScreenVisualize::Apply() {
    if (!m_window && !m_window_disabled) {
        CreateGlfwWindow(Name());
    }
    if (m_window) {
        // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        // do all memcpy
        cudaStreamSynchronize(m_cuda_stream);
        if (m_bufferR8) {
            cudaMemcpyAsync(m_hostR8->Buffer.get(), m_bufferR8->Buffer.get(),
                            m_bufferR8->Width * m_bufferR8->Height * sizeof(char), cudaMemcpyDeviceToHost,
                            m_cuda_stream);
        } else if (m_bufferRGBA8) {
            cudaMemcpyAsync(m_hostRGBA8->Buffer.get(), m_bufferRGBA8->Buffer.get(),
                            m_hostRGBA8->Width * m_hostRGBA8->Height * sizeof(PixelRGBA8), cudaMemcpyDeviceToHost,
                            m_cuda_stream);
        } else if (m_bufferDI) {
            cudaMemcpyAsync(m_hostDI->Buffer.get(), m_bufferDI->Buffer.get(),
                            m_hostDI->Width * m_hostDI->Height * sizeof(PixelDI), cudaMemcpyDeviceToHost,
                            m_cuda_stream);
        } else if (m_bufferRangeRcs) {
            cudaMemcpyAsync(m_hostRangeRcs->Buffer.get(), m_bufferRangeRcs->Buffer.get(),
                            m_hostRangeRcs->Width * m_hostRangeRcs->Height * sizeof(PixelRangeRcs),
                            cudaMemcpyDeviceToHost, m_cuda_stream);
        } else {
            throw std::runtime_error("No buffer incoming for visualization");
        }

        // lock the glfw mutex because from here on out, we don't want to be interrupted
        std::lock_guard<std::mutex> lck(s_glfwMutex);

        // do window prep stuff
        glfwMakeContextCurrent(m_window.get());
        if (!m_gl_tex_id) {
            glGenTextures(1, &m_gl_tex_id);
            glBindTexture(GL_TEXTURE_2D, m_gl_tex_id);

            // Change these to GL_LINEAR for super- or sub-sampling
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

            // GL_CLAMP_TO_EDGE for linear filtering, not relevant for nearest.
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }
        glBindTexture(GL_TEXTURE_2D, m_gl_tex_id);

        // Set Viewport to window dimensions
        int window_w, window_h;
        glfwGetWindowSize(m_window.get(), &window_w, &window_h);
        glViewport(0, 0, window_w, window_h);

        // update the textures, making sure data has finished memcpy first
        cudaStreamSynchronize(m_cuda_stream);
        if (m_bufferR8) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, m_hostR8->Width, m_hostR8->Height, 0, GL_RED, GL_UNSIGNED_BYTE,
                         m_hostR8->Buffer.get());
        } else if (m_bufferRGBA8) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_hostRGBA8->Width, m_hostRGBA8->Height, 0, GL_RGBA,
                         GL_UNSIGNED_BYTE, m_hostRGBA8->Buffer.get());
        } else if (m_bufferDI) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, m_hostDI->Width, m_hostDI->Height, 0, GL_RG, GL_FLOAT,
                         m_hostDI->Buffer.get());
            // TODO: support dual returns
        } else if (m_bufferRangeRcs) {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, m_hostRangeRcs->Width, m_hostRangeRcs->Height, 0, GL_RG, GL_FLOAT,
                         m_hostRangeRcs->Buffer.get());
            // TODO: support dual returns
        } else {
            throw std::runtime_error("No buffer incoming for visualization");
        }

        // update the window

        // 1:1 texel to pixel mapping with glOrtho(0, 1, 0, 1, -1, 1) setup:
        // The quad coordinates go from lower left corner of the lower left pixel
        // to the upper right corner of the upper right pixel.
        // Same for the texel coordinates.

        glEnable(GL_TEXTURE_2D);
        glBegin(GL_QUADS);
        glTexCoord2f(0.0f, 0.0f);
        glVertex2f(0.0f, 0.0f);
        glTexCoord2f(1.0f, 0.0f);
        glVertex2f(1.0f, 0.0f);
        glTexCoord2f(1.0f, 1.0f);
        glVertex2f(1.0f, 1.0f);
        glTexCoord2f(0.0f, 1.0f);
        glVertex2f(0.0f, 1.0f);
        glEnd();
        glDisable(GL_TEXTURE_2D);

        glfwSwapBuffers(m_window.get());
        glfwPollEvents();
        // std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> wall_time = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        // std::cout << "Render to window time: " << wall_time.count() << std::endl;
    }
}
CH_SENSOR_API void ChFilterFullScreenVisualize::Initialize(std::shared_ptr<ChSensor> pSensor,
                                                 std::shared_ptr<SensorBuffer>& bufferInOut) {
    if (!bufferInOut)
        InvalidFilterGraphNullBuffer(pSensor);

    auto pOptixSen = std::dynamic_pointer_cast<ChOptixSensor>(pSensor);
    if (!pOptixSen) {
        InvalidFilterGraphSensorTypeMismatch(pSensor);
    }
    m_cuda_stream = pOptixSen->GetCudaStream();
    m_bufferR8 = std::dynamic_pointer_cast<SensorDeviceR8Buffer>(bufferInOut);
    m_bufferRGBA8 = std::dynamic_pointer_cast<SensorDeviceRGBA8Buffer>(bufferInOut);
    m_bufferDI = std::dynamic_pointer_cast<SensorDeviceDIBuffer>(bufferInOut);
    m_bufferRangeRcs = std::dynamic_pointer_cast<SensorDeviceRangeRcsBuffer>(bufferInOut);

    if (m_bufferR8) {
        m_hostR8 = chrono_types::make_shared<SensorHostR8Buffer>();
        std::shared_ptr<char[]> b(cudaHostMallocHelper<char>(m_bufferR8->Width * m_bufferR8->Height),
                                  cudaHostFreeHelper<char>);
        m_hostR8->Buffer = std::move(b);
        m_hostR8->Width = m_bufferR8->Width;
        m_hostR8->Height = m_bufferR8->Height;

    } else if (m_bufferRGBA8) {
        m_hostRGBA8 = chrono_types::make_shared<SensorHostRGBA8Buffer>();
        std::shared_ptr<PixelRGBA8[]> b(cudaHostMallocHelper<PixelRGBA8>(m_bufferRGBA8->Width * m_bufferRGBA8->Height),
                                        cudaHostFreeHelper<PixelRGBA8>);
        m_hostRGBA8->Buffer = std::move(b);
        m_hostRGBA8->Width = m_bufferRGBA8->Width;
        m_hostRGBA8->Height = m_bufferRGBA8->Height;
    } else if (m_bufferDI) {
        m_hostDI = chrono_types::make_shared<SensorHostDIBuffer>();
        std::shared_ptr<PixelDI[]> b(cudaHostMallocHelper<PixelDI>(m_bufferDI->Width * m_bufferDI->Height),
                                     cudaHostFreeHelper<PixelDI>);
        m_hostDI->Buffer = std::move(b);
        m_hostDI->Width = m_bufferDI->Width;
        m_hostDI->Height = m_bufferDI->Height;
    } else if (m_bufferRangeRcs) {
        m_hostRangeRcs = chrono_types::make_shared<SensorHostRangeRcsBuffer>();
        std::shared_ptr<PixelRangeRcs[]> b(
            cudaHostMallocHelper<PixelRangeRcs>(m_bufferRangeRcs->Width * m_bufferRangeRcs->Height),
            cudaHostFreeHelper<PixelRangeRcs>);
        m_hostRangeRcs->Buffer = std::move(b);
        m_hostRangeRcs->Width = m_bufferRangeRcs->Width;
        m_hostRangeRcs->Height = m_bufferRangeRcs->Height;
    } else {
        InvalidFilterGraphBufferTypeMismatch(pSensor);
    }
}

CH_SENSOR_API void ChFilterFullScreenVisualize::CreateGlfwWindow(std::string window_name) {
    // if we've already made the window, there's nothing to do.
    if (m_window)
        return;

    OnNewWindow();  // OnNewWindow will need to lock inside itself

    std::lock_guard<std::mutex> lck(s_glfwMutex);
		if (m_is_fullscreen)
		    m_window.reset(
						glfwCreateWindow(static_cast<GLsizei>(m_w), static_cast<GLsizei>(m_h), window_name.c_str(), glfwGetPrimaryMonitor(), NULL));
		else
				m_window.reset(
						glfwCreateWindow(static_cast<GLsizei>(m_w), static_cast<GLsizei>(m_h), window_name.c_str(), NULL, NULL));
    if (m_window) {
        glfwMakeContextCurrent(m_window.get());
        glfwSwapInterval(0);  // disable vsync as we are "fast as possible"
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, 1, 0, 1, -1, 1);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        glViewport(0, 0, m_w, m_h);
    } else {
        std::cerr << "WARNING: requested window could not be created by GLFW. Will proceed with no window.\n";
        m_window_disabled = true;
    }
}

CH_SENSOR_API void ChFilterFullScreenVisualize::OnNewWindow() {
    std::lock_guard<std::mutex> l(s_glfwMutex);
    if (s_windowCount++ == 0) {
        glfwInit();
    }
}
CH_SENSOR_API void ChFilterFullScreenVisualize::OnCloseWindow() {
    std::lock_guard<std::mutex> l(s_glfwMutex);
    if (--s_windowCount == 0) {
        glfwTerminate();
    }
}

}  // namespace sensor
}  // namespace chrono
