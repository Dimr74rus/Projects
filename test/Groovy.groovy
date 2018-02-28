@Grab(group = 'cc.redberry', module = 'groovy', version = '1.1.11-SNAPSHOT')
import cc.redberry.groovy.Redberry

import cc.redberry.core.context.CC
import cc.redberry.core.context.OutputFormat
import cc.redberry.core.tensor.Power
import cc.redberry.core.tensor.Product
import cc.redberry.core.tensor.Tensor
import cc.redberry.core.transformations.Transformation
import cc.redberry.groovy.Redberry

import static cc.redberry.core.utils.TensorUtils.Exponent
import static cc.redberry.core.utils.TensorUtils.info
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

    def $Degree = { Tensor expr -> Exponent('a'.t, expr) }
    def $Select = { Tensor expr, exponent ->
        (expr.class == Product || expr.class == Power) && ($Degree(expr) > exponent) ? 0 : expr
    }
    def Select = { exponent ->
        def innerTransform = { Tensor expr -> expr.transformParentAfterChild { subExpr -> $Select(subExpr, exponent) } } as Transformation
        { Tensor expr -> (ExpandAndEliminate['d^\\alpha_\\alpha = 4'.t & innerTransform] & innerTransform) >> expr } as Transformation
    }

    // превый порядок(т.к. по умолчанию в пакете g - метрика, то g сдесь соответствует \eta в pdf файле с вычислениями, и наоборот)
    def n = 'n_{\\alpha \\beta} = g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]'.t
    def Dg = 'Dg_{\\alpha\\beta} = D[x^\\gamma][g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]] * a * \\ksi^{\\gamma}[x^\\alpha]'.t
    def n1 = 'n1_{\\alpha \\beta} = n_{ \\alpha \\beta} + Dg_{ \\alpha \\beta}'.t
    def dx1 = 'dx1^{\\alpha}_{\\beta} =  D[x^\\beta][x^{\\alpha} + a * \\ksi^{\\alpha}[x^\\alpha]]'.t
    def n2 = 'n2_{\\alpha \\beta} = dx1^{\\gamma}_{\\alpha} * dx1^{\\sigma}_{\\beta} * n1_{\\gamma \\sigma}'.t

    n1 <<= Dg & n

    n2 <<= dx1 & n1 & Dg & n
    n2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t

    def res1 = 'res1^{\\gamma}_{\\mu \\nu}'.t.eq((n2 & Differentiate['x_{\\gamma}']) >> 'n2_{\\mu \\nu}'.t)
    res1 = 'res1_{\\beta}'.t.eq(res1 >> 'res1_{\\mu \\beta}^{\\mu}'.t)
    res1 <<= EliminateMetrics & 'd^\\alpha_\\alpha = 4'.t
    println n1
    println res1

    def tt = 'a**2'.t
    tt.parentAfterChild { a -> print a.toString() + ', ' }

    res1 = Select(1) >> res1

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
    //n2 = Select(4) >> n2

    def res2 = 'res2_{\\mu\\nu\\gamma}'.t.eq((n2 & Differentiate['x^{\\gamma}']) >> 'n2_{\\mu \\nu}'.t)
    res2 = 'res2_{\\beta}'.t.eq((res2 & dx2) >> 'dx2^{\\gamma \\alpha} * res2_{ \\alpha \\beta \\gamma}'.t)
    res2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t
    //res2 = Select(2) >> res2
    res2 = 'res2_{\\alpha}'.t.eq((res2 & res1) >> 'res2_{\\alpha} - res1_{\\alpha}'.t)

    println res2
    println res2.toString(OutputFormat.UTF8)
}
