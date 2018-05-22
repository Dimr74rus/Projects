import cc.redberry.groovy.Redberry

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.Power
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry
import cc.redberry.core.indices.IndicesSymmetries

import static cc.redberry.core.utils.TensorUtils.Exponent
//import static cc.redberry.core.utils.TensorUtils.info
import static cc.redberry.groovy.RedberryStatic.*

use(Redberry) {
    def I = 'I_{\\nu \\rho} = (q_{\\nu} n_{\\alpha} n_{\\beta} + n_{\\nu} k_{\\alpha} n_{\\beta} + n_{\\nu} n_{\\alpha} (q_\\beta - k_\\beta))*(g^{\\alpha \\delta} - (1-s) k^\\alpha k^\\delta/(k^\\gamma k^\\gamma))*(k_\\delta n_\\rho n_\\sigma + n_\\delta q_\\rho n_\\sigma + n_\\delta n_\\rho (q_\\sigma - k_\\sigma))*(g^{\\sigma \\beta} - (1-s) (q^{\\alpha} - k^{\\alpha}) (q^{\\beta} - k^\\delta)/((q^\\gamma-k^\\gamma) (q^\\gamma - k^\\gamma)))'.t
    I << EliminateMetrics & 'd^\\alpha_\\alpha = 4'.t

    println I
    //println res2.toString(OutputFormat.UTF8)
}
