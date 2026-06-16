module CirculatorySystemModels

using ModelingToolkit
using Setfield

export Pin, OnePort, Ground, Resistor, QResistor, PoiseuilleResistor, Capacitor, Inductance, Compliance, Elastance, VariableElastance, ConstantPressure, ConstantFlow, DrivenPressure, DrivenFlow, DHChamber, ShiChamber, ShiAtrium, ShiHeart, WK3, WK3E, CR, CRL, RRCR, ShiSystemicLoop, ShiPulmonaryLoop, ResistorDiode, OrificeValve, ShiValve, MynardValve_SemiLunar, MynardValve_Atrioventricular


@independent_variables t

D = Differential(t)

function Pin(; name)
        @variables p(t) q(t) [connect = Flow]
        sys = System(Equation[], t, [p, q], []; name)
        sys = @set sys.connector_type = ModelingToolkit.connector_type(sys)
        return sys
end


function Ground(; name, P=0.0)
        @named g = Pin()
        ps = @parameters P = P
        eqs = [g.p ~ P]
        compose(System(eqs, t, [], ps; name), g)
end


function OnePort(; name)
        @named out = Pin()
        @named in = Pin()
        sts = @variables begin
                Δp(t)
                q(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                0 ~ in.q + out.q
                q ~ in.q
        ]
        compose(System(eqs, t, sts, []; name), in, out)
end

function OnePortWithExtPressure(; name)
        @named out = Pin()
        @named in = Pin()
        @named ep = Pin()
        sts = @variables begin
                Δp(t)
                q(t)
                pg(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                0 ~ in.q + out.q
                0 ~ ep.q
                q ~ in.q
                pg ~ p - ep.p
        ]
        compose(System(eqs, t, sts, []; name), in, out, ep)
end


"""
`Resistor(;name, R=1.0)`

Implements the resistor using Ohm's law to represent a vessels linear resistance to blood flow.

Parameter is in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg
`q` calculated in cm^3/s (ml/s)

Named parameters:

`R`:       Resistance of the vessel to the fluid in mmHg*s/ml
"""
function Resistor(; name, R=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters R = R
        eqs = [Δp ~ -q * R]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`QResistor(;name, K=1.0)`

Implements the quadratic resistor to represent a vessels non-linear resistance to blood flow.

Parameters are in the cm, g, s system.
Pressures in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`K`: non-linear resistance of the vessel to the fluid in mmHg*s^2/ml^2
"""
function QResistor(; name, K=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters K = K
        eqs = [Δp ~ -q * abs(q) * K]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`Capacitor(;name, C=1.0)`

Implements a capacitor to represent vessel capacitance.

Parameters are in the cm, g, s system.
Pressures in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`C`:      capacitance of the vessel in ml/mmHg
"""
function Capacitor(; name, C=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters C = C
        eqs = [D(Δp) ~ -q / C]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`Inductance(;name, L=1.0)`

Implements the inductance to represent blood inertance.

Parameters are in the cm, g, s system.
Pressures in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`L`:       Inertia of the fluid in mmHg*s^2/ml
"""
function Inductance(; name, L=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters L = L
        eqs = [D(q) ~ -Δp / L]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`PoiseuilleResistor(;name, μ=3e-2, r=0.5, L=5)`

Implements the resistance following the Poiseuille law.

Parameters are in the cm, g, s system.
Pressures in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`μ`:       viscosity of fluid in dyne s / cm^2

`r`:       radius of vessel segmenty in cm

`L`:       length of vessel segment in cm
"""
function PoiseuilleResistor(; name, μ=3e-2, r=0.1, L=1)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters μ = μ r = r L = L
        R = 8 * μ * L / (π * r^4) * (1 / 1333.2)
        eqs = [Δp ~ -q * R]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`Compliance(; name, V₀=0.0, C=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)`

Implements the compliance of a vessel.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`p` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`V₀`:               Unstressed volume ml

`C`:                Vessel compliance in ml/mmHg


`inP`:             (Bool) formulate in dp/dt (default: false)

`has_ep`:          (Bool) if true, add a parameter `p₀` for pressure offset
                   e.g., for thoracic pressure (default: false)

`p₀`:              External pressure in mmHg (e.g., thorax pressure, default: 0.0)
                   _Note: if this argument is set, it will be used, even if `has_ep` is
                   `false`. `has_ep` only controls if `p₀` will be exposed as a parameter!_

has_variable_ep`: (Bool) expose pin for variable external pressure (default: false)
                   This pin can be connected to another pin or function providing external pressure.
                   _Note: if `has_variable_ep` is set to `true` this pin is created, independent of
                   `has_ep`!_
"""
function Compliance(; name, V₀=0.0, C=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)
        @named in = Pin()
        @named out = Pin()

        if has_variable_ep
                @named ep = Pin()
        end

        sts = @variables begin
                V(t)
                p(t)
        end

        ps = @parameters begin
                V₀ = V₀
                C = C
        end

        # Add the thoracic pressure variant

        D = Differential(t)

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if has_variable_ep
                push!(sts,
                        (@variables p_rel(t))[1]
                )
                if has_ep
                        push!(ps,
                                (@parameters p₀ = p₀)[1]
                        )
                end
                push!(eqs,
                        p_rel ~ ep.p + p₀,
                        ep.q ~ 0
                )
        elseif has_ep
                push!(ps,
                        (@parameters p₀ = p₀)[1]
                )
                p_rel = p₀
        else
                p_rel = p₀
        end

        if inP
                push!(eqs,
                        V ~ (p - p_rel) * C + V₀,
                        D(p) ~ (in.q + out.q) * 1 / C
                )
        else
                push!(eqs,
                        p ~ (V - V₀) / C + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        if has_variable_ep
                compose(System(eqs, t, sts, ps; name=name), in, out, ep)
        else
                compose(System(eqs, t, sts, ps; name=name), in, out)
        end
end


"""
`Elastance(; name, V₀=0.0, E=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)`

Implements the elastance of a vessel. Elastance more commonly used to describe the heart.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`p` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`V₀`:             Unstressed volume ml

`E`:              Vessel elastance in ml/mmHg. Equivalent to compliance as E=1/C

`inP`:            (Bool) formulate in dp/dt (default: false)

`has_ep`:         (Bool) if true, add a parameter `p₀` for pressure offset
                  e.g., for thoracic pressure (default: false)

`p₀`:             External pressure in mmHg (e.g., thorax pressure, default: 0.0)
                  _Note: if this argument is set, it will be used, even if `has_ep` is
                  `false`. `has_ep` only controls if `p₀` will be exposed as a parameter!_

has_variable_ep`: (Bool) expose pin for variable external pressure (default: false)
                   This pin can be connected to another pin or function providing external pressure.
                   _Note: if `has_variable_ep` is set to `true` this pin is created, independent of
                   `has_ep`!_
"""
function Elastance(; name, V₀=0.0, E=1.0, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)
        @named in = Pin()
        @named out = Pin()

        if has_variable_ep
                @named ep = Pin()
        end

        sts = @variables begin
                V(t)
                p(t)
        end

        ps = @parameters begin
                V₀ = V₀
                E = E
        end

        D = Differential(t)

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if has_variable_ep
                push!(sts,
                        (@variables p_rel(t))[1]
                )
                if has_ep
                        push!(ps,
                                (@parameters p₀ = p₀)[1]
                        )
                end
                push!(eqs,
                        p_rel ~ ep.p + p₀,
                        ep.q ~ 0
                )
        elseif has_ep
                push!(ps,
                        (@parameters p₀ = p₀)[1]
                )
                p_rel = p₀
        else
                p_rel = p₀
        end

        if inP
                push!(eqs,
                        V ~ (p - p_rel) / E + V₀,
                        D(p) ~ (in.q + out.q) * E
                )
        else
                push!(eqs,
                        p ~ (V - V₀) * E + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        if has_variable_ep
                compose(System(eqs, t, sts, ps; name=name), in, out, ep)
        else
                compose(System(eqs, t, sts, ps; name=name), in, out)
        end
end


"""
`VariableElastance(; name, V₀=0.0, C=1.0, Escale=1.0, fun, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)`

`VariableElastance` is defined based on the `Elastance` element,
but has a time varying elastance function modelling
the contraction of muscle fibres.

Named parameters:

`V₀`:              stress-free volume (zero pressure volume)

`Escale`:          scaling factor (elastance factor)

`fun`:             function object for elastance (must be `fun(t)`)

`inP`:             (Bool) formulate in dp/dt (default: false)

`has_ep`:          (Bool) if true, add a parameter `p₀` for pressure offset
                   e.g., for thoracic pressure (default: false)

`p₀`:              External pressure in mmHg (e.g., thorax pressure, default: 0.0)
                   _Note: if this argument is set, it will be used, even if `has_ep` is
                   `false`. `has_ep` only controls if `p₀` will be exposed as a parameter!_

has_variable_ep`: (Bool) expose pin for variable external pressure (default: false)
                   This pin can be connected to another pin or function providing external pressure.
                   _Note: if `has_variable_ep` is set to `true` this pin is created, independent of
                   `has_ep`!_
"""
function VariableElastance(; name, V₀=0.0, C=1.0, Escale=1.0, fun, inP=false, has_ep=false, has_variable_ep=false, p₀=0.0)
        @named in = Pin()
        @named out = Pin()

        if has_variable_ep
                @named ep = Pin()
        end

        sts = @variables begin
                V(t)
                p(t)
        end

        ps = @parameters begin
                V₀ = V₀
                C = C
        end

        D = Differential(t)
        E = Escale * fun(t)

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if has_variable_ep
                push!(sts,
                        (@variables p_rel(t))[1]
                )
                push!(eqs,
                        p_rel ~ ep.p,
                        ep.q ~ 0
                )
        elseif has_ep
                push!(ps,
                        (@parameters p₀ = p₀)[1]
                )
                p_rel = p₀
        else
                p_rel = p₀
        end

        if inP
                push!(eqs,
                        V ~ (p - p_rel) / E + V₀,
                        D(p) ~ (in.q + out.q) * E + V * D(E(t))
                )
        else
                push!(eqs,
                        p ~ (V - V₀) * E + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        if has_variable_ep
                compose(System(eqs, t, sts, ps; name=name), in, out, ep)
        else
                compose(System(eqs, t, sts, ps; name=name), in, out)
        end
end


"""
`ConstantPressure(;name, P=1.0)`

Implements a constant pressure source to a system.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`P`:     Constant pressure in mmHg
"""
function ConstantPressure(; name, P=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters P = P
        eqs = [Δp ~ P]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`ConstantFlow(;name, Q=1.0)`

Implements a constant flow source to a system.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`Q`:     Constant flow in cm^3/s (ml/s)
"""
function ConstantFlow(; name, Q=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters Q = Q
        eqs = [q ~ Q]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`DrivenPressure(; name, fun, P=1.0)`

Implements a driven pressure source to a system modulated by a function provided.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`P`:     Constant pressure in mmHg

`fun`:   Function which modulates the input
"""
function DrivenPressure(; name, fun=sin, P=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters P = P
        eqs = [Δp ~ P * fun(t)]
        extend(System(eqs, t, [], ps; name), oneport)
end

"""
`DrivenFlow(; name, fun, Q=1.0)`

Implements a driven flow source to a system.

Parameters are in the cm, g, s system.
Pressure in mmHg.
`Δp` is calculated in mmHg,
`q` is calculated in cm^3/s (ml/s).

Named parameters:

`Q`:     Constant flow in cm^3/s (ml/s).

`τ`     Length of cardiac cycle is s

`fun`:   Function which modulates the input
"""
function DrivenFlow(; name, fun=sin, Q=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters Q = Q
        eqs = [q ~ Q * fun(t)]
        extend(System(eqs, t, [], ps; name), oneport)
end

"""
`DHChamber(;name, V₀, Eₘᵢₙ, n₁, n₂, τ, τ₁, τ₂, k, Eshift=0.0, inP=false)`

The Double Hill chamber/ventricle model is defined based on the vessel
element, but has a time varying elastance function modelling the contraction
of muscle fibres

The time varying elastance is calculated using the Double Hill model.

This model uses external helper functions `elastance` and `delastance`
which describe the elastance function and the first derivative of it.

It calculates the elastance as:

E(t) = (Eₘₐₓ - Eₘᵢₙ) * e(t) + Eₘᵢₙ

where e(t) is the Double-Hill function.

Named parameters:

`V₀`:     stress-free volume (zero pressure volume)

`p₀`      pressure offset (defaults to zero)
          this is present in some papers (e.g. Shi), so is
          provided here for conformity. Defaults to 0.0

`Eₘᵢₙ`:   minimum elastance

`Eₘₐₓ`:   maximum elastance

`n₁`:     rise coefficient

`n₂`:     fall coefficient

`τ`:      pulse length [s]

`τ₁`:     rise timing parameter[s]

`τ₂`:     fall timimg paramter [s]

`k`:      elastance factor*

`Eshift`: time shift of contraction (for atria)

`inP`:    (Bool) formulate in dp/dt (default: false)

*Note: `k` is not an independent parameter, it is a scaling factor that corresponds
to 1/max(e(t)), which ensures that e(t) varies between zero and 1.0, such that
E(t) varies between Eₘᵢₙ and Eₘₐₓ.
"""
function DHChamber(; name, V₀, p₀=0.0, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, k, Eshift=0.0, inP=false)
        @named in = Pin()
        @named out = Pin()
        sts = @variables begin
                V(t)
                p(t)
        end
        ps = @parameters begin
                V₀ = V₀
                p₀ = p₀
                Eₘᵢₙ = Eₘᵢₙ
                Eₘₐₓ = Eₘₐₓ
                n₁ = n₁
                n₂ = n₂
                τ = τ
                τ₁ = τ₁
                τ₂ = τ₂
                k = k
                Eshift = Eshift
        end

        E = DHelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
        DE = DHdelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
        p_rel = p₀

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if inP
                push!(eqs,
                        V ~ (p - p_rel) / E + V₀,
                        D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE
                )
        else
                push!(eqs,
                        p ~ (V - V₀) * E + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        compose(System(eqs, t, sts, ps; name), in, out)
end


"""
`DHelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)`

Helper function for `DHChamber`
"""
function DHelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
        tᵢ = rem(t + (1 - Eshift) * τ, τ)

        return (Eₘₐₓ - Eₘᵢₙ) * k * ((tᵢ / τ₁)^n₁ / (1 + (tᵢ / τ₁)^n₁)) * (1 / (1 + (tᵢ / τ₂)^n₂)) + Eₘᵢₙ
end


"""
`DHdelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)`

Helper function for `DHChamber`
"""
function DHdelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
        tᵢ = rem(t + (1 - Eshift) * τ, τ)

        de = ((Eₘₐₓ - Eₘᵢₙ) * k * (τ₂^n₂ * n₁ * tᵢ^(n₁ - 1) *
                                   (τ₁^n₁ * τ₂^n₂ + τ₁^n₁ * tᵢ^n₂ + τ₂^n₂ * tᵢ^n₁ + tᵢ^(n₁ + n₂)) -
                                   τ₂^n₂ * tᵢ^n₁ * (τ₁^n₁ * n₂ * tᵢ^(n₂ - 1) + τ₂^n₂ * n₁ * tᵢ^(n₁ - 1) +
                                                    (n₁ + n₂) * tᵢ^(n₁ + n₂ - 1))) /
              (τ₁^n₁ * τ₂^n₂ + τ₁^n₁ * tᵢ^n₂ + τ₂^n₂ * tᵢ^n₁ + tᵢ^(n₁ + n₂))^2)

        return de
end


"""
`ShiChamber(;name, V₀, p₀=0.0, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift=0.0)`

Implemention of a ventricle following Shi/Korakianitis.

This model uses external helper function `shiElastance`
which describes the elastance function.

Named parameters:

`V₀`     stress-free volume (zero pressure volume)

`p₀`     pressure offset (defaults to zero)
         this is present in the original paper, so is
         provided here for conformity. Defaults to 0.0

`Eₘᵢₙ`   minimum elastance

`τ`      pulse length

`τₑₛ`    end systolic time (end of rising cosine)

`τₑₚ`    end pulse time (end of falling cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle

`inP`:    (Bool) formulate in dp/dt (default: false)
"""
function ShiChamber(; name, V₀, p₀=0.0, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift, inP=false)
        @named in = Pin()
        @named out = Pin()
        sts = @variables begin
                V(t)
                p(t)
        end
        ps = @parameters begin
                V₀ = V₀
                p₀ = p₀
                Eₘᵢₙ = Eₘᵢₙ
                Eₘₐₓ = Eₘₐₓ
                τ = τ
                τₑₛ = τₑₛ
                τₑₚ = τₑₚ
                Eshift = Eshift
        end

        E = ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
        DE = DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)
        p_rel = p₀

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if inP
                push!(eqs,
                        V ~ (p - p_rel) / E + V₀,
                        D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE
                )
        else
                push!(eqs,
                        p ~ (V - V₀) * E + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        compose(System(eqs, t, sts, ps; name), in, out)
end


"""
`ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)`

Elastance function `E(t)` for ventricle simulation based on Shi's
double cosine function.

Parameters:

`Eₘᵢₙ`: minimum elastance (diastole)

`Eₘₐₓ`: maximum elastance (systole)

`τₑₛ`: end systolic time (end of rising cosine)

`τₑₚ`: end of pulse time (end of falling cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle
"""
function ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

        tᵢ = rem(t + (1 - Eshift) * τ, τ)

        Eₚ = (tᵢ <= τₑₛ) * (1 - cos(tᵢ / τₑₛ * pi)) / 2 +
             (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * (1 + cos((tᵢ - τₑₛ) / (τₑₚ - τₑₛ) * pi)) / 2 +
             (tᵢ <= τₑₚ) * 0

        E = Eₘᵢₙ + (Eₘₐₓ - Eₘᵢₙ) * Eₚ

        return E
end


"""
DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

Helper function for `ShiChamber`

Derivative of the elastance function `E(t)` for ventricle simulation based on Shi's
double cosine function.

Parameters:

`Eₘᵢₙ`: minimum elastance (diastole)

`Eₘₐₓ`: maximum elastance (systole)

`τₑₛ`: end systolic time (end of rising cosine)

`τₑₚ`: end of pulse time (end of falling cosine)

`Eshift`: time shift of contraction (for atria), set to `0` for ventricle
"""
function DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, τₑₛ, τₑₚ, Eshift)

        tᵢ = rem(t + (1 - Eshift) * τ, τ)

        DEₚ = (tᵢ <= τₑₛ) * pi / τₑₛ * sin(tᵢ / τₑₛ * pi) / 2 +
              (tᵢ > τₑₛ) * (tᵢ <= τₑₚ) * pi / (τₑₚ - τₑₛ) * sin((τₑₛ - tᵢ) / (τₑₚ - τₑₛ) * pi) / 2
        (tᵢ <= τₑₚ) * 0
        DE = (Eₘₐₓ - Eₘᵢₙ) * DEₚ

        return DE
end


"""
`ShiAtrium(;name, V₀, p₀, Eₘᵢₙ, Eₘₐₓ, τ, τpwb, τpww, inP = false)`

Implementation of the Atrium following Shi/Korakianitis.

Named parameters:

name    name of the element

`V₀`    Unstressed chamber volume in ml

`p₀`    Unstressed chamber pressure in mmHg

`Eₘᵢₙ`  Minimum elastance (diastole) in mmHg/ml

`Eₘₐₓ`  Maximum elastance (systole) in mmHg/ml

`τ`     Length of cardiac cycle in s

`τpwb`  Atrial contraction time in s

`τpww`  Atrial offset time in s
"""
function ShiAtrium(; name, V₀, p₀, Eₘᵢₙ, Eₘₐₓ, τ, τpwb, τpww, inP=false)
        @named in = Pin()
        @named out = Pin()
        sts = @variables begin
                V(t)
                p(t)
        end
        ps = @parameters begin
                V₀ = V₀
                p₀ = p₀
                Eₘᵢₙ = Eₘᵢₙ
                Eₘₐₓ = Eₘₐₓ
                τ = τ
                τpwb = τpwb
                τpww = τpww
        end

        # adjust timing parameters to fit the elastance functions for the ventricle
        # define elastance based on ventricle E function
        E = ShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, 0.5 * τpww, τpww, τpwb)
        DE = DShiElastance(t, Eₘᵢₙ, Eₘₐₓ, τ, 0.5 * τpww, τpww, τpwb)
        p_rel = p₀

        eqs = [
                0 ~ in.p - out.p
                p ~ in.p
        ]

        if inP
                push!(eqs,
                        V ~ (p - p_rel) / E + V₀,
                        D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE
                )
        else
                push!(eqs,
                        p ~ (V - V₀) * E + p_rel,
                        D(V) ~ in.q + out.q
                )
        end

        compose(System(eqs, t, sts, ps; name), in, out)
end


"""
`ShiHeart(; name, τ, LV.V₀, LV.p0, LV.Emin, LV.Emax, LV.τes, LV.τed, LV.Eshift,
        RV.V₀, RV.p0, RV.Eₘᵢₙ, RV.Eₘₐₓ, RV.τes, RV.τed, RV.Eshift, LA.V₀, LA_p0,
        LA.Emin, LA.Emax, LA.τes, LA.τed, LA.Eshift, RA.V₀, RA.p0, RA.Emin,
        RA.Emax, RA.τes, RA.τed, RA.Eshift, AV.CQ, AV.Kp, AV.Kf, AV.Kb, AV.Kv,
        AV.θmax, AV.θmin, PV.CQ, PV.Kp, PV.Kf, PV.Kb, PV.Kv, PV.θmax, PV.θmin,
        MV.CQ, MV.Kp, MV.Kf, MV.Kb, MV.Kv, MV.θmax, MV.θmin, TV.CQ, TV.Kp, TV.Kf,
        TV.Kb, TV.Kv, TV.θmax, TV.θmin)`

Models a whole heart, made up of 2 ventricles (Left & Right Ventricle) and 2 atria (Left & Right atrium)
created from the ShiChamber element. Includes the 4 corresponding valves (Aortic, Mitral, Pulmonary and Tricuspid valve) created using the ShiValve element.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).
Maximum and Minimum angles given in rad, to convert from degrees multiply angle by pi/180.

Named parameters:

`τ`         Length of the cardiac cycle in s

`LV_V₀`     Unstressed left ventricular volume in ml

`LV_p0`     Unstressed left ventricular pressure in mmHg

`LV_Emin`   Minimum left ventricular elastance (diastole) in mmHg/ml

`LV_Emax`   Maximum left ventricular elastance (systole) in mmHg/ml

`LV_τes`    Left ventricular end systolic time in s

`LV_τed`    Left ventricular end distolic time in s

`LV_Eshift` Shift time of contraction - 0 for left ventricle

`RV_V₀`     Unstressed right ventricular volume in ml

`RV_p0`     Unstressed right ventricular pressure in mmHg

`RV_Emin`   Minimum right ventricular elastance (diastole) in mmHg/ml

`RV_Emax`   Maximum right ventricular elastance (systole) in mmHg/ml

`RV_τes`    Right ventricular end systolic time in s

`RV_τed`    Right ventricular end distolic time in s

`RV_Eshift` Shift time of contraction - 0 for right ventricle

`LA_V₀`     Unstressed left atrial volume in ml

`LA_p0`     Unstressed left atrial pressure in mmHg

`LA_Emin`   Minimum left atrial elastance (diastole) in mmHg/ml

`LA_Emax`   Maximum left atrial elastance (systole) in mmHg/ml

`LA_τes`    Left atrial end systolic time in s

`LA_τed`    Left atrial end distolic time in s

`LA_Eshift` Shift time of contraction in s

`RA_V₀`     Unstressed right atrial volume in ml

`RA_p0`     Unstressed right atrial pressure in mmHg

`RA_Emin`   Minimum right atrial elastance (diastole) in mmHg/ml

`RA_Emax`   Maximum right atrial elastance (systole) in mmHg/ml

`RA_τes`    Right atrial end systolic time in s

`RA_τed`    Right atrial end distolic time in s

`RA_Eshift` Shift time of contraction in s

`AV_CQ`     Aortic valve flow coefficent in ml/(s*mmHg^0.5)

`AV_Kp`     Pressure effect on the aortic valve in rad/(s^2*mmHg)

`AV_Kf`     Frictional effect on the aortic valve in 1/s

`AV_Kb`     Fluid velocity effect on the aortic valve in rad/(s*m)

`AV_Kv`     Vortex effect on the aortic valve in rad/(s*m)

`AV_θmax`   Aortic valve maximum opening angle in rad

`AV_θmin`   Aortic valve minimum opening angle in rad

`MV_CQ`     Mitral valve flow coefficent in ml/(s*mmHg^0.5)

`MV_Kp`     Pressure effect on the mitral valve in rad/(s^2*mmHg)

`MV_Kf`     Frictional effect on the mitral valve in 1/s

`MV_Kb`     Fluid velocity effect on the mitral valve in rad/(s*m)

`MV_Kv`     Vortex effect on the mitral valve in rad/(s*m)

`MV_θmax`   Mitral valve maximum opening angle in rad

`MV_θmin`   Mitral valve minimum opening angle in rad

`PV_CQ`     Pulmonary valve flow coefficent in ml/(s*mmHg^0.5)

`PV_Kp`     Pressure effect on the pulmonary valve in rad/(s^2*mmHg)

`PV_Kf`     Frictional effect on the pulmonary valve in 1/s

`PV_Kb`     Fluid velocity effect on the pulmonary valve in rad/(s*m)

`PV_Kv`     Vortex effect on the pulmonary valve in rad/(s*m)

`PV_θmax`   Pulmonary valve maximum opening angle in rad

`PV_θmin`   Pulmonary valve minimum opening angle in rad

`TV_CQ`     Tricuspid valve flow coefficent in ml/(s*mmHg^0.5)

`TV_Kp`     Pressure effect on the tricuspid valve in rad/(s^2*mmHg)

`TV_Kf`     Frictional effect on the tricuspid valve in 1/s

`TV_Kb`     Fluid velocity effect on the tricuspid valve in rad/(s*m)

`TV_Kv`     Vortex effect on the pulmonary valve in rad/(s*m)

`TV_θmax`   Tricuspid valve maximum opening angle in rad

`TV_θmin`   Tricuspid valve minimum opening angle in rad
"""
function ShiHeart(; name, τ,
                LV_V₀, LV_p₀=0.0, LV_Eₘᵢₙ, LV_Eₘₐₓ, LV_τ, LV_τₑₛ, LV_τₑₚ, LV_Eshift=0.0,
                RV_V₀, RV_p₀=0.0, RV_Eₘᵢₙ, RV_Eₘₐₓ, RV_τ, RV_τₑₛ, RV_τₑₚ, RV_Eshift=0.0,
                LA_V₀, LA_p₀=0.0, LA_Eₘᵢₙ, LA_Eₘₐₓ, LA_τ, LA_τₑₛ, LA_τₑₚ, LA_Eshift=0.0,
                RA_V₀, RA_p₀=0.0, RA_Eₘᵢₙ, RA_Eₘₐₓ, RA_τ, RA_τₑₛ, RA_τₑₚ, RA_Eshift=0.0,
                AV_CQ, AV_Kp, AV_Kf, AV_Kb, AV_Kv, AV_θmax, AV_θmin,
                MV_CQ, MV_Kp, MV_Kf, MV_Kb, MV_Kv, MV_θmax, MV_θmin,
                TV_CQ, TV_Kp, TV_Kf, TV_Kb, TV_Kv, TV_θmax, TV_θmin,
                PV_CQ, PV_Kp, PV_Kf, PV_Kb, PV_Kv, PV_θmax, PV_θmin)
        @named LHin = Pin()
        @named LHout = Pin()
        @named RHin = Pin()
        @named RHout = Pin()
        @named in = Pin()
        @named out = Pin()
        sts = @variables begin
                Δp(t)
                q(t)
        end

        # Ventricles and atria
        @named LV = ShiChamber(V₀=LV_V₀, p₀=LV_p₀, Eₘᵢₙ=LV_Eₘᵢₙ, Eₘₐₓ=LV_Eₘₐₓ, τ=LV_τ, τₑₛ=LV_τₑₛ, τₑₚ=LV_τₑₚ, Eshift=LV_Eshift)
        @named RV = ShiChamber(V₀=RV_V₀, p₀=RV_p₀, Eₘᵢₙ=RV_Eₘᵢₙ, Eₘₐₓ=RV_Eₘₐₓ, τ=RV_τ, τₑₛ=RV_τₑₛ, τₑₚ=RV_τₑₚ, Eshift=RV_Eshift)
        @named LA = ShiChamber(V₀=LA_V₀, p₀=LA_p₀, Eₘᵢₙ=LA_Eₘᵢₙ, Eₘₐₓ=LA_Eₘₐₓ, τ=LA_τ, τₑₛ=LA_τₑₛ, τₑₚ=LA_τₑₚ, Eshift=LA_Eshift)
        @named RA = ShiChamber(V₀=RA_V₀, p₀=RA_p₀, Eₘᵢₙ=RA_Eₘᵢₙ, Eₘₐₓ=RA_Eₘₐₓ, τ=RA_τ, τₑₛ=RA_τₑₛ, τₑₚ=RA_τₑₚ, Eshift=RA_Eshift)
        # Valves
        @named AV = ShiValve(CQ=AV_CQ, Kp=AV_Kp, Kf=AV_Kf, Kb=AV_Kb, Kv=AV_Kv, θmax=AV_θmax, θmin=AV_θmin)
        @named MV = ShiValve(CQ=MV_CQ, Kp=MV_Kp, Kf=MV_Kf, Kb=MV_Kb, Kv=MV_Kv, θmax=MV_θmax, θmin=MV_θmin)
        @named TV = ShiValve(CQ=TV_CQ, Kp=TV_Kp, Kf=TV_Kf, Kb=TV_Kb, Kv=TV_Kv, θmax=TV_θmax, θmin=TV_θmin)
        @named PV = ShiValve(CQ=PV_CQ, Kp=PV_Kp, Kf=PV_Kf, Kb=PV_Kb, Kv=PV_Kv, θmax=PV_θmax, θmin=PV_θmin)

        eqs = [
                Δp ~ out.p - in.p
                q ~ in.q
                connect(LHin, LA.in)
                connect(LA.out, MV.in)
                connect(MV.out, LV.in)
                connect(LV.out, AV.in)
                connect(AV.out, LHout)
                connect(RHin, RA.in)
                connect(RA.out, TV.in)
                connect(TV.out, RV.in)
                connect(RV.out, PV.in)
                connect(PV.out, RHout)
        ]

        compose(System(eqs, t, sts, []; name), LHin, LHout, RHin, RHout, in, out, LV, RV, LA, RA, AV, MV, TV, PV)
end


"""
`ResistorDiode(;name, R=1e-3)`

Implements the resistance across a valve following Ohm's law exhibiting diode like behaviour.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)

Named parameters:

`R`     Resistance across the valve in mmHg*s/ml
"""
function ResistorDiode(; name, R=1e-3)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters R = R
        eqs = [q ~ -Δp / R * (Δp < 0)]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`OrificeValve(;name, CQ=1.0)`

Implements the square-root pressure-flow relationship across a valve.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)

Named parameters:

`CQ`    Flow coefficent in ml/(s*mmHg^0.5)
"""
function OrificeValve(; name, CQ=1.0)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters CQ = CQ
        eqs = [q ~ (Δp < 0) * CQ * sqrt(sign(Δp) * Δp)]
        extend(System(eqs, t, [], ps; name), oneport)
end


"""
`ShiValve(; name, CQ, Kp, Kf, Kb, Kv, θmax, θmin)`

Implements the Shi description for valve opening and closing, full description in [Shi].

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)
Maximum and Minimum angles given in rad, to convert from degrees multiply angle by pi/180.

Named parameters:

`CQ`    Flow coefficent in ml/(s*mmHg^0.5)

`Kp`    Pressure effect on the valve in rad/(s^2*mmHg)

`Kf`    Frictional effect on the valve in 1/s

`Kb`    Fluid velocity effect on the valve in rad/(s*m)

`Kv`    Vortex effect on the valve in rad/(s*m)

`θmax`  Valve maximum opening angle in rad

`θmin`  Valve minimum opening angle in rad
"""
function ShiValve(; name, CQ, Kp, Kf, Kb, Kv, θmax, θmin)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters CQ = CQ Kp = Kp Kf = Kf Kb = Kb Kv = Kv θmax = θmax θmin = θmin
        sts = @variables θ(t) ω(t) AR(t) Fp(t) Ff(t) Fb(t) Fv(t) F(t)
        D = Differential(t)

        function bb_affect_min!(mod, obs, integ, ctx)
                return (; ω = 0.0)
        end
        function bb_affect_max!(mod, obs, integ, ctx)
                return (; ω = 0.0)
        end
        limits = [
                [θ ~ θmax] => (bb_affect_max!, (; ω))
                [θ ~ θmin] => (bb_affect_min!, (; ω))
        ]

        # make θmax the real opening angle and define a θmaxopen for a healthy valve
        # that means we can use θmax as a stenosis parameter
        θmaxopen = 75 * pi / 180

        eqs = [
                # Forces/Moments
                Fp ~ Kp * -Δp * cos(θ)                 # pressure
                Ff ~ -Kf * ω                           # friction
                Fb ~ Kb * q * cos(θ)                       # Fluid Velocity
                Fv ~ -Kv * q * (q > 0) * sin(2θ)       # vortex behind leaflets
                F ~ Fp + Ff + Fb + Fv                  # total force/moment on leaflets
                #ODEs
                D(θ) ~ ω
                D(ω) ~ F * ((θ < θmax) * (F > 0) + (θ > θmin) * (F < 0))
                # Opening ratio
                #AR ~ ((1 - cos(θ))^2) / ((1 - cos(θmax))^2)
                AR ~ ((1 - cos(θ))^2) / ((1 - cos(θmaxopen))^2)
                # Flow equation
                q ~ -sign(Δp) * CQ * AR * sqrt(abs(Δp))
        ]

        # include the `continuous_events` definition `limits` in the ODE system
        # this is the MTK equivalent to callbacks
        extend(System(eqs, t, sts, ps; name=name, continuous_events=limits), oneport)
end


"""
`MynardValve_SemiLunar(; name, ρ, Leff, Mrg, Mst, Ann, Kvc, Kvo)`

Implements the Mynard description for flow across the semilunar valves, full description in [Mynard].
This valve description corresponds to the semilunar valves where interia is an effect we consider.

Note: The minimum level of regurgitation has to be set to machine precision eps()

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)
p is scaled to ensure units are consistent throughout.

Named parameters:

name    name of the element
`ρ`     Blood density in g/cm^3
`Leff`  An effective length in cm
`Mrg`   Level of regurgitation exhibited by a valve in DN
`Mst`   Level of stenosis exhibited by a valve in DN
`Ann`   Annulus area in cm^2
`Kvc`   Valve closing rate coefficent in  cm^2/(dynes*s)
`Kvo`   Valve opening rate coefficent in cm^2/(dynes*s)

p is calculated in mmHg
q is calculated in cm^3/s (ml/s)
"""
function MynardValve_SemiLunar(; name, ρ, Leff, Mrg, Mst, Ann, Kvc, Kvo)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters ρ = ρ Leff = Leff Mrg = Mrg Mst = Mst Ann = Ann Kvc = Kvc Kvo = Kvo
        sts = @variables begin
                Aeff(t)
                ζ(t)
                B(t)
                Aeff_min(t)
                Aeff_max(t)
                L(t)
        end
        Δp = -1333.22 * Δp
        eqs = [
                # Opening ratio
                D(ζ) ~ (Δp > 0) * ((1 - ζ) * Kvo * Δp) + (Δp < 0) * (ζ * Kvc * Δp)
                Aeff_min ~ Mrg * Ann + eps()
                Aeff_max ~ Mst * Ann
                Aeff ~ (Aeff_max - Aeff_min) * ζ + Aeff_min
                # Flow equation
                B ~ ρ / (2 * Aeff^2)
                L ~ ρ * Leff / Aeff
                D(q) ~ (Δp - B * q * abs(q)) * 1 / L
        ]
        extend(System(eqs, t, sts, ps; name), oneport)
end


"""
`MynardValve_Atrioventricular(; name, ρ, Mrg, Mst, Ann, Kvc, Kvo)`

Implements the Mynard description for flow across the atrioventricular valves, full description in [Mynard].
This valve description corresponds to the atrioventricular valves where interia is not considered.

Note: The minimum level of regurgitation has to be set to machine precision eps()

Parameters are in the cm, g, s system.
Pressure in mmHg.
Flow in cm^3/s (ml/s)
p is scaled to ensure units are consistent throughout.

Named parameters:

name    name of the element
`ρ`     Blood density in g/cm^3
`Mrg`   Level of regurgitation exhibited by a valve in DN
`Mst`   Level of stenosis exhibited by a valve in DN
`Ann`   Annulus area in cm^2
`Kvc`   Valve closing rate coefficent in  cm^2/(dynes*s)
`Kvo`   Valve opening rate coefficent in cm^2/(dynes*s)

p is calculated in mmHg
q is calculated in cm^3/s (ml/s)
"""
function MynardValve_Atrioventricular(; name, ρ, Mrg, Mst, Ann, Kvc, Kvo)
        @named oneport = OnePort()
        @unpack Δp, q = oneport
        ps = @parameters ρ = ρ Mrg = Mrg Mst = Mst Ann = Ann Kvc = Kvc Kvo = Kvo
        sts = @variables begin
                Aeff(t)
                ζ(t)
                B(t)
                Aeff_min(t)
                Aeff_max(t)
                L(t)
        end
        p = -1333.22 * p
        eqs = [
                # Opening ratio
                D(ζ) ~ (Δp > 0) * ((1 - ζ) * Kvo * Δp) + (Δp < 0) * (ζ * Kvc * Δp)
                Aeff_min ~ Mrg * Ann + eps()
                Aeff_max ~ Mst * Ann
                Aeff ~ (Aeff_max - Aeff_min) * ζ + Aeff_min
                # Flow equation
                B ~ ρ / (2 * Aeff^2)
                q ~ sqrt(1 / B * abs(Δp)) * sign(Δp)
        ]
        extend(System(eqs, t, sts, ps; name), oneport)
end


"""
`WK3(;name, Rc=1.0, Rp=1.0, C=1.0)`

Implements the 3 element windkessel model.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`Rp`:      Peripheral resistance in mmHg*s/ml

`C`:       Arterial compliance in ml/mmHg
"""
function WK3(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named C = Capacitor(C=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, Rc.in)
                connect(Rc.out, Rp.in, C.in)
                connect(Rp.out, C.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, C, ground), oneport)
end


"""
`WK3E(;name, Rc=1.0, Rp=1.0, E=1.0)`

Implements the 3 element windkessel model. With a vessel elastance instead of a capacitor.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`Rp`:      Peripheral resistance in mmHg*s/ml

`E`:       Arterial elastance in mmHg/ml
"""
function WK3E(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named E = Elastance(E=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, Rc.in)
                connect(Rc.out, E.in)
                connect(E.out, Rp.in)
                connect(Rp.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, E, ground), oneport)
end


"""
`WK4_S(;name, Rc=1.0, L=1.0, Rp=1.0, C=1.0)`

Implements the 4 element windkessel model with serial inertance.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`L`:       Inertance/Inductance in mmHg*s^2*ml^-1

`Rp`:      Peripheral resistance in mmHg*s/ml

`C`:       Arterial compliance in ml/mmHg
"""
function WK4_S(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named C = Capacitor(C=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, Rc.in)
                connect(Rc.out, L.in)
                connect(L.out, C.in, Rp.in)
                connect(Rp.out, C.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, C, L, ground), oneport)
end


"""
`WK4_SE(;name, Rc=1.0, L=1.0, Rp=1.0, E=1.0)`

Implements the 4 element windkessel model with serial inertance.
With a vessel elastance instead of a capacitor.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`L`:       Inertance/Inductance in mmHg*s^2*ml^-1

`Rp`:      Peripheral resistance in mmHg*s/ml

`E`:       Arterial elastance in mmHg/ml
"""
function WK4_SE(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named E = Elastance(E=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, Rc.in)
                connect(Rc.out, L.in)
                connect(L.out, E.in)
                connect(E.out, Rp.in)
                connect(Rp.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, E, L, ground), oneport)
end


"""
`WK4_P(;name, Rc=1.0, L=1.0, Rp=1.0, C=1.0)`

Implements the 4 element windkessel model with parallel inertance.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`L`:       Inertance/Inductance in mmHg*s^2*ml^-1

`Rp`:      Peripheral resistance in mmHg*s/ml

`C`:       Arterial compliance in ml/mmHg
"""
function WK4_P(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named C = Capacitor(C=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, L.in, Rc.in)
                connect(L.out, Rc.out, C.in, Rp.in)
                connect(Rp.out, C.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, C, L, ground), oneport)
end


"""
`WK4_PE(;name, Rc=1.0, L=1.0, Rp=1.0, E=1.0)`

Implements the 4 element windkessel model with parallel inertance.
With a vessel elastance instead of a capacitor.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s)

Named parameters:

`Rc`:      Characteristic impedance in mmHg*s/ml

`L`:       Inertance/Inductance in mmHg*s^2*ml^-1

`Rp`:      Peripheral resistance in mmHg*s/ml

`E`:       Arterial elastance in mmHg/ml
"""
function WK4_PE(; name)
        @named oneport = OnePort()
        @named Rc = Resistor(R=1.0)
        @named Rp = Resistor(R=1.0)
        @named E = Elastance(E=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, L.in, Rc.in)
                connect(L.out, Rc.out, E.in)
                connect(E.out, Rp.in)
                connect(Rp.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), Rc, Rp, E, L, ground), oneport)
end


function WK5(; name)
        @named oneport = OnePort()
        @named R1 = Resistor(R=1.0)
        @named C1 = Capacitor(C=1.0)
        @named R2 = Resistor(R=1.0)
        @named C2 = Capacitor(C=1.0)
        @named R3 = Resistor(R=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, R1.in)
                connect(R1.out, C1.in, R2.in)
                connect(R2.out, C2.in, R3.in)
                connect(R3.out, C1.out, C2.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), R1, C1, R2, C2, R3, L, ground), oneport)
end


function WK5E(; name)
        @named oneport = OnePort()
        @named R1 = Resistor(R=1.0)
        @named E1 = Elastance(E=1.0)
        @named R2 = Resistor(R=1.0)
        @named E2 = Elastance(E=1.0)
        @named R3 = Resistor(R=1.0)
        @named L = Inductance(L=1.0)
        @named ground = Ground()
        eqs = [
                connect(oneport.in, R1.in)
                connect(R1.out, E1.in)
                connect(E1.out, R2.in)
                connect(R2.out, E2.in)
                connect(E2.out, R3.in)
                connect(R3.out, oneport.out)
        ]
        extend(compose(System(eqs, t, [], []; name), R1, E1, R2, E2, R3, L, ground), oneport)
end


"""
`CR(;name, R=1.0, C=1.0)`

Implements the compliance, resistor subsystem.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`R`:       Component resistance in mmHg*s/ml

`C`:       Component compliance in ml/mmHg
"""
function CR(; name, R=1.0, C=1.0)
        @named in = Pin()
        @named out = Pin()
        @named R = Resistor(R=R)
        @named C = Compliance(C=C)
        sts = @variables begin
                Δp(t)
                q(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                q ~ in.q
                connect(in, C.in)
                connect(C.out, R.in)
                connect(R.out, out)
        ]
        compose(System(eqs, t, sts, []; name), in, out, R, C)
end


"""
`CRL(;name, C=1.0, R=1.0, L=1.0)`

Implements the compliance, resistor, inductance subsystem.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`C`:       Component compliance in ml/mmHg

`R`:       Component resistance in mmHg*s/ml

`L`:       Component blood inertia in mmHg*s^2/ml
"""
function CRL(; name, C=1.0, R=1.0, L=1.0)
        @named in = Pin()
        @named out = Pin()
        @named C = Compliance(C=C)
        @named R = Resistor(R=R)
        @named L = Inductance(L=L)
        sts = @variables begin
                Δp(t)
                q(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                q ~ in.q
                connect(in, C.in)
                connect(C.out, R.in)
                connect(R.out, L.in)
                connect(L.out, out)
        ]
        compose(System(eqs, t, sts, []; name), in, out, C, R, L)
end


"""
`RRCR(;name, R1=1.0, R2=1.0, R3=1.0, C=1.0)`

Implements the resistor, resistor, compliance, resistor subsystem.

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`R1`:      Component resistance in mmHg*s/ml

`R2`:      Component resistance in mmHg*s/ml

`C`:       Component compliance in ml/mmHg

`R3`:      Component resistance in mmHg*s/ml
"""
function RRCR(; name)
        @named in = Pin()
        @named out = Pin()
        @named ep = Pin()
        @named R1 = Resistor(R=1.0)
        @named R2 = Resistor(R=1.0)
        @named C = Compliance(C=1.0)
        @named R3 = Resistor(R=1.0)
        sts = @variables begin
                p(t)
                q(t)
        end
        eqs = [
                p ~ out.p - in.p
                q ~ in.q
                connect(in, R1.in)
                connect(R1.out, R2.in)
                connect(R2.out, C.in)
                connect(C.out, R3.in)
                connect(R3.out, out)
                connect(C.ep, ep)
        ]
        compose(System(eqs, t, sts, []; name), in, out, ep, R1, R2, C, R3)
end


"""
`ShiSystemicLoop(; name, SAS.C, SAS.R, SAS.L, SAT.C, SAT.R, SAT.L, SAR.R, SCP.R, SVN.C, SVN.R)`

Implements systemic loop as written by Shi in [Shi].

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`SAS_C`:   Aortic sinus compliance in ml/mmHg

`SAS_R`:   Aortic sinus resistance in mmHg*s/ml

`SAS_L`:   Aortic sinus inductance in mmHg*s^2/ml

`SAT_C`:   Artery compliance in ml/mmHg

`SAT_R`:   Artery resistance in mmHg*s/ml

`SAT_L`:   Artery inductance in mmHg*s^2/ml

`SAR_R`:   Arteriole resistance in mmHg*s/ml

`SCP_R`:   Capillary resistance in mmHg*s/ml

`SVN_C`:   Vein compliance in ml/mmHg

`SVN_R`:   Vein resistance in mmHg*s/ml
"""
function ShiSystemicLoop(; name,
                SAS_C=1.0, SAS_R=1.0, SAS_L=1.0,
                SAT_C=1.0, SAT_R=1.0, SAT_L=1.0,
                SAR_R=1.0, SCP_R=1.0,
                SVN_C=1.0, SVN_R=1.0)
        @named in = Pin()
        @named out = Pin()
        # These are the components the subsystem is made of:
        ## Systemic Aortic Sinus ##
        @named SAS = CRL(C=SAS_C, R=SAS_R, L=SAS_L)
        ## Systemic Artery ##
        @named SAT = CRL(C=SAT_C, R=SAT_R, L=SAT_L)
        ## Systemic Arteriole ##
        @named SAR = Resistor(R=SAR_R)
        ## Systemic Capillary ##
        @named SCP = Resistor(R=SCP_R)
        ## Systemic Vein ##
        @named SVN = CR(C=SVN_C, R=SVN_R)
        sts = @variables begin
                Δp(t)
                q(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                q ~ in.q
                connect(in, SAS.in)
                connect(SAS.out, SAT.in)
                connect(SAT.out, SAR.in)
                connect(SAR.out, SCP.in)
                connect(SCP.out, SVN.in)
                connect(SVN.out, out)
        ]
        compose(System(eqs, t, sts, []; name), in, out, SAS, SAT, SAR, SCP, SVN)
end


"""
`ShiPulmonaryLoop(; name, PAS.C, PAS.R, PAS.L, PAT.C, PAT.R, PAT.L, PAR.R, PCP.R, PVN.C, PVN.R)`

Implements systemic loop as written by Shi in [Shi].

Parameters are in the cm, g, s system.
Pressure in mmHg.
Volume in ml.
Flow in cm^3/s (ml/s).

Named parameters:

`PAS__C`:   Artery sinus compliance in ml/mmHg

`PAS__R`:   Artery sinus resistance in mmHg*s/ml

`PAS__L`:   Artery sinus Inductance in mmHg*s^2/ml

`PAT__C`:   Artery compliance in ml/mmHg

`PAT__R`:   Artery resistance in mmHg*s/ml

`PAT__L`:   Artery Inductance in mmHg*s^2/ml

`PAR__R`:   Arteriole resistance in mmHg*s/ml

`PCP__R`:   Capillary resistance in mmHg*s/ml

`PVN__C`:   Vein compliance in ml/mmHg

`PVN__R`:   Vein resistance in mmHg*s/ml
"""
function ShiPulmonaryLoop(; name,
                PAS_C=1.0, PAS_R=1.0, PAS_L=1.0,
                PAT_C=1.0, PAT_R=1.0, PAT_L=1.0,
                PAR_R=1.0, PCP_R=1.0,
                PVN_C=1.0, PVN_R=1.0)
        @named in = Pin()
        @named out = Pin()
        # These are the components the subsystem is made of:
        ## Pulmonary Aortic Sinus ##
        @named PAS = CRL(C=PAS_C, R=PAS_R, L=PAS_L)
        ## Pulmonary Artery ##
        @named PAT = CRL(C=PAT_C, R=PAT_R, L=PAT_L)
        ## Pulmonary Arteriole ##
        @named PAR = Resistor(R=PAR_R)
        ## Pulmonary Capillary ##
        @named PCP = Resistor(R=PCP_R)
        ## Pulmonary Vein ##
        @named PVN = CR(C=PVN_C, R=PVN_R)
        sts = @variables begin
                Δp(t)
                q(t)
        end
        eqs = [
                Δp ~ out.p - in.p
                q ~ in.q
                connect(in, PAS.in)
                connect(PAS.out, PAT.in)
                connect(PAT.out, PAR.in)
                connect(PAR.out, PCP.in)
                connect(PCP.out, PVN.in)
                connect(PVN.out, out)
        ]
        compose(System(eqs, t, sts, []; name), in, out, PAS, PAT, PAR, PCP, PVN)
end


end
