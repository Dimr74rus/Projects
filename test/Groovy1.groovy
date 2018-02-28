/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2018:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.groovy.scripts

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

// подключены сразу все модули, возможно не все они нужны
use(Redberry) {
    // функция, которая отбрасывает порядки выше od у выражения expr
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

    def $Degree = { expr -> Exponent(expr, 'a'.t) }
    def $Select = { expr, exponent ->
        (expr.class == Product || expr.class == Power) && $Degree(expr) > exponent ? 0.t : expr
    }
    def Select = { exponent ->
        def innerTransform = { Tensor expr -> expr.transformParentAfterChild { subExpr -> $Select(subExpr, exponent) } } as Transformation
        { Tensor expr -> (ExpandAndEliminate['d^\\alpha_\\alpha = 4'.t & innerTransform] & innerTransform) >> expr } as Transformation
    }


    //CC.parserAllowsSameVariance = true

    // превый порядок(т.к. по умолчанию в пакете g - метрика, то g сдесь соответствует \eta в pdf файле с вычислениями, и наоборот)
    def n = 'n_{\\alpha \\beta}  = g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]'.t
    def Dg = 'Dg_{\\alpha\\beta} = D[x^\\gamma][g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\alpha]] * a * \\ksi^{\\gamma}[x^\\alpha]'.t
    def n1 = 'n1_{\\alpha \\beta} = n_{ \\alpha \\beta} + Dg_{ \\alpha \\beta}'.t
    def dx1 = 'dx1^{\\alpha}_{\\beta} =  D[x^\\beta][x^{\\alpha} + a * \\ksi^{\\alpha}[x^\\alpha]]'.t
    def n2 = 'n2_{\\alpha \\beta} = dx1^{\\gamma}_{\\alpha} * dx1^{\\sigma}_{\\beta} * n1_{\\gamma \\sigma}'.t

    n2 <<= dx1 & n1 & Dg & n
    n2 <<= dx1 & n1 & Dg & n
    n2 <<= ExpandAndEliminate & 'd^\\alpha_\\alpha = 4'.t

    println info(n2)


    def res1 = 'res1_\\gamma\\mu\\nu'.t.eq((n2 & Differentiate['x^{\\gamma}']) >> 'n2_{\\mu \\nu}'.t)
    res1 = res1 >> 'res1_{\\mu \\beta \\mu}'.t
    res1 <<= EliminateMetrics & 'd^\\alpha_\\alpha = 4'.t
    println n2
    println res1

    res1 = Select(1) >> res1

    println res1.toString(OutputFormat.LaTeX)

    // второй порядок
    n = 'g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\delta]'.t
    Dg = 'D[x^\\gamma][g_{\\alpha \\beta} + a * h_{\\alpha \\beta}[x^\\psi]] * (a * \\ksi^{\\gamma}[x^\\delta] + a**2 * \\sigma^{\\gamma}[x^\\delta])'.t
    DDg = "D[x^\\epsilon][($Dg)] * (a * \\ksi^{\\epsilon}[x^\\delta] + a**2 * \\sigma^{\\epsilon}[x^\\delta])".t
    n1 = 'n_{ \\alpha \\beta} + Dg_{ \\alpha \\beta} + DDg_{ \\alpha \\beta}'.t
    dx1 = 'D[x^\\beta][x^{\\alpha} + a * \\ksi^{\\alpha}[x^\\alpha] + a**2 * \\sigma^{\\alpha}[x^\\alpha]]'.t
    def dx2 = 'd^{\\gamma}_{\\delta} - a * D[x^\\delta][\\ksi^{\\gamma}[x^\\alpha]] + a**2 * (D[x^\\psi][\\ksi^{\\gamma}[x^\\alpha]] * D[x^\\delta][\\ksi^{\\psi}[x^\\alpha]] - D[x^\\delta][\\sigma^{\\gamma}[x^\\alpha]])'.t
    n2 = 'dx1^{\\gamma}_{\\alpha} * dx1^{\\sigma}_{\\beta} * n1_{\\gamma \\sigma}'.t

    n1 = ('Dg_{\\alpha \\beta}'.t.eq(Dg) & 'DDg_{\\alpha \\beta}'.t.eq(DDg) & 'n_{ \\alpha \\beta}'.t.eq(n)) >> n1

    n2 = ('n1_{ \\alpha \\beta}'.t.eq(n1) & 'dx1^{\\alpha}_{\\beta}'.t.eq(dx1) & 'n_{ \\alpha \\beta}'.t.eq(n)) >> n2
    println 'before select'
    n2 = Select(4) >> n2
    println 'after select'

    def res2 = Differentiate['x^{\\gamma}'] >> n2
    res2 = ('res_{ \\alpha \\beta \\gamma}'.t.eq(res2) & 'dx2^{\\gamma}_{\\delta}'.t.eq(dx2)) >> 'g^{\\delta \\alpha} * (dx2^{\\gamma}_{\\delta}) * res_{ \\alpha \\beta \\gamma}'.t
    res2 = Select(2) >> res2
    res2 = "$res2 - $res1".t
    res2 <<= ExpandAndEliminate
    println info(res2)

    println res2.toString(OutputFormat.UTF8)
}