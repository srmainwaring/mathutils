//
// Created by frongere on 15/12/17.
//

#ifndef FRYDOM_INTEGRATION1D_H
#define FRYDOM_INTEGRATION1D_H

#include <vector>

namespace mathutils {

    template <class Scalar>
    class Integrate1d {

        enum INTEGRATION_DATA {
            VECTOR,
            FUNCTION
        };

    protected:
        Scalar m_Xmin = 0.;
        Scalar m_Xmax = 1.;
        unsigned int m_NbPoints;

        std::vector<Scalar> m_Y;
        Scalar (*m_integrand)(Scalar x) = nullptr; // Function to be integrated

        INTEGRATION_DATA m_DataSource= VECTOR;

        bool c_Computed = false;

        Scalar c_Result = 0.;


    private:
        virtual void Compute() = 0;

    public:
        Integrate1d(Scalar (*F)(Scalar x), Scalar xmin, Scalar xmax, unsigned int nbPoints)
                : m_Xmin(xmin),
                  m_Xmax(xmax),
                  m_NbPoints(nbPoints),
                  c_Computed(false) {

            m_integrand = F;
            m_DataSource = FUNCTION;
        }

        Integrate1d(const std::vector<Scalar>& Y, Scalar xmin, Scalar xmax, unsigned int nbPoints) :
                m_Xmin(xmin),
                m_Xmax(xmax),
                m_NbPoints(nbPoints),
                c_Computed(false) {  // FIXME: pas besoin de mettre nbPoints...

            m_Y = Y;
            m_DataSource = VECTOR;
        }

        Scalar Get() {
            if (!c_Computed) {
                Compute();
                c_Computed = true;
            }
            return c_Result;
        }

    };

    template <class Scalar>
    class Trapz : public Integrate1d {

    private:

        void Compute() {
            // Implementation of the trapezoidal method for 1d numerical integration

            Scalar dx = (m_Xmax-m_Xmin) / Scalar(m_NbPoints);

            if (m_DataSource == FUNCTION) {
                m_Y.clear();
                m_Y.swap(m_Y);

                Scalar x(m_Xmin);

                for (unsigned int i=0; i<m_NbPoints; i++) {
                    m_Y.push_back((*m_integrand)(x));
                    x += dx;
                }

            }

            Scalar sum = 0.0;
            Scalar x(m_Xmin + dx);

            for (unsigned int i=1; i<m_NbPoints; i++) {
                sum += m_Y[i];
                x += dx;
            }

            c_Result = m_Y[0] + m_Y[m_Y.size()] + 2.0*sum;
            c_Result *= 0.5 * dx;

        }

    public:
        Trapz(Scalar (*F)(Scalar x), Scalar xmin, Scalar xmax, unsigned int nbPoints)
                : Integrate1d(F, xmin, xmax, nbPoints) {}

        Trapz(const std::vector<Scalar>& Y, Scalar xmin, Scalar xmax, unsigned int nbPoints)
                : Integrate1d(Y, xmin, xmax, nbPoints) {}

    };









}  // end namespace mathutils


#endif //FRYDOM_INTEGRATION1D_H
