@Grab(group = 'cc.redberry', module = 'groovy', version = '1.1.11-SNAPSHOT')
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
    /*
    def trc = { Tensor expr, od ->
        for (int i = 15; i > od; --i) {
            def t = expr
            def k = 1
            for (int j = i; j > 0; j--) {
                t <<= Differentiate['a']
                k *= j
            }
            expr -= 'a'.t**i * t / k
        }
        expr <<= Expand & CollectScalars
        expr
    }
    */
    def $Degree = { expr -> Exponent(expr, 'a'.t) }
    def $Select = { expr, exponent ->
        (expr.class == Product || expr.class == Power) && $Degree(expr) > exponent ? 0.t : expr
    }
    def Select = { exponent ->
        def innerTransform = { Tensor expr -> expr.transformParentAfterChild { subExpr -> $Select(subExpr, exponent) } } as Transformation
        { Tensor expr -> (ExpandAndEliminate['d^\\alpha_\\alpha = 4'.t & innerTransform] & innerTransform) >> expr } as Transformation
    }

    //addSymmetry('h_{\\alpha \\beta}'.t, [1, 0].p)
    //'n_{\\alpha \\beta}'.t.indices.symmetries.addSymmetry( [1, 0].p )

    // превый порядок(т.к. по умолчанию в пакете g - метрика, то g сдесь соответствует \eta в pdf файле с вычислениями, и наоборот)
    def n = 'n_{\\alpha \\beta} = g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]'.t
    def Dg = 'Dg_{\\alpha\\beta} = D[x^\\gamma][g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]] * a * \\ksi^{\\gamma}[x^\\alpha]'.t
    def n1 = 'n1_{\\alpha \\beta} = n_{ \\alpha \\beta} + Dg_{ \\alpha \\beta}'.t
    def dx1 = 'dx1^{\\alpha}_{\\beta} =  D[x^\\beta][x^{\\alpha} + a * \\ksi^{\\alpha}[x^\\alpha]]'.t
    def n2 = 'n2_{\\alpha \\beta} = dx1^{\\gamma}_{\\alpha} * dx1^{\\sigma}_{\\beta} * n1_{\\gamma \\sigma}'.t

    def symmetries = findIndicesSymmetries('_{\\alpha \\beta}'.si, 'h_{\\alpha \\beta}'.t)
    for (s in symmetries)
        println s

    n1 <<= Dg & n

    def t = 'b*c'.t
    println Exponent(t, 'a'.t)

    n2 <<= dx1 & n1 & Dg & n
    n2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t

    def res1 = 'res1^{\\gamma}_{\\mu \\nu}'.t.eq((n2 & Differentiate['x_{\\gamma}']) >> 'n2_{\\mu \\nu}'.t)
    res1 = 'res1_{\\beta}'.t.eq(res1 >> 'res1_{\\mu \\beta}^{\\mu}'.t)
    res1 <<= EliminateMetrics & 'd^\\alpha_\\alpha = 4'.t

    //def tt = 'a**8*T_{\\alpha}'.t
    //tt.parentAfterChild { a -> print a.toString() + ', ' }

    res1 = Select(1) >> res1
    println res1

    // второй порядок
    n = 'n_{\\alpha \\beta} = g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\delta]'.t
    def Vx1 = 'Vx1^{\\gamma}[x^\\delta] = a * \\ksi^{\\gamma}[x^\\delta] + a**2 * \\sigma^{\\gamma}[x^\\delta]'.t
    Dg = 'Dg_{\\alpha\\beta}[x^\\gamma] = D[x^\\gamma][g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\psi]] * Vx1^{\\gamma}[x^\\delta]'.t
    DDg = 'DDg_{\\alpha\\beta} =  D[x^\\epsilon][Dg_{\\alpha\\beta}[x^\\gamma]] * Vx1^{\\epsilon}[x^\\delta]'.t
    n1 = 'n1_{ \\alpha \\beta} = n_{ \\alpha \\beta} + Dg_{ \\alpha \\beta}[x^\\gamma] + DDg_{ \\alpha \\beta}'.t
    dx1 = 'dx1^{\\alpha}_{\\beta} = D[x^\\beta][x^{\\alpha} + a * \\ksi^{\\alpha}[x^\\alpha] + a**2 * \\sigma^{\\alpha}[x^\\alpha]]'.t
    def dx2 = 'dx2^{\\gamma}_{\\delta} = d^{\\gamma}_{\\delta} - a * D[x^\\delta][\\ksi^{\\gamma}[x^\\alpha]] + a**2 * (D[x^\\psi][\\ksi^{\\gamma}[x^\\alpha]] * D[x^\\delta][\\ksi^{\\psi}[x^\\alpha]] - D[x^\\delta][\\sigma^{\\gamma}[x^\\alpha]])'.t
    n2 = 'n2_{\\alpha \\beta} = dx1^{\\gamma}_{\\alpha} * dx1^{\\sigma}_{\\beta} * n1_{\\gamma \\sigma}'.t

    n1 <<= DDg & Dg & n & Vx1
    println n1
    println n1.toString(OutputFormat.UTF8)

    n2 <<= n1 & dx1 & n & Vx1
    n2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t
    n2 = Select(4) >> n2

    def res2 = 'res2_{\\mu\\nu\\gamma}'.t.eq((n2 & Differentiate['x^{\\gamma}']) >> 'n2_{\\mu \\nu}'.t)
    res2 = 'res2_{\\beta}'.t.eq((res2 & dx2) >> 'dx2^{\\gamma \\alpha} * res2_{ \\alpha \\beta \\gamma}'.t)
    res2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t
    res2 = Select(1) >> res2
    res2 = 'res2^{\\alpha}'.t.eq((res2 & res1) >> 'res2^{\\alpha} - res1^{\\alpha}'.t)

    println res2

    //def I = 'I_{\\nu \\rho} = (q_{\\nu} n_{\\alpha} n_{\\beta} + n_{\\nu} k_{\\alpha} n_{\\beta} + n_{\\nu} n_{\\alpha} (q_\\beta - k_\\beta))*(g^{\\alpha \\delta} - (1-s) k^\\alpha k^\\delta/(k^\\gamma k^\\gamma))*(k_\\delta n_\\rho n_\\sigma + n_\\delta q_\\rho n_\\sigma + n_\\delta n_\\rho (q_\\sigma - k_\\sigma))*(g^{\\sigma \\beta} - (1-s) (q^{\\alpha} - k^{\\alpha}) (q^{\\beta} - k^\\delta)/((q^\\gamma-k^\\gamma) (q^\\gamma - k^\\gamma)))'.t
    //I << EliminateMetrics & 'd^\\alpha_\\alpha = 4'.t

    //println I
    //println res2.toString(OutputFormat.UTF8)
}
