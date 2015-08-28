using Distributions

immutable NormalInverseGaussian <: ContinuousUnivariateDistribution
  μ::Float64
  α::Float64
  β::Float64
  δ::Float64

  function NormalInverseGaussian(μ::Real, α::Real, β::Real, δ::Real)
    new(μ, α, β, δ)
  end
end

params(d::NormalInverseGaussian) = (d.μ, d.α, d.β, d.δ)

function logpdf(d::NormalInverseGaussian, x::Float64)
  μ, α, β, δ = params(d)
  log(α*δ) + log(besselK(α*sqrt(δ^2+(x-μ)^2),1)) - log(π*sqrt(δ^2+(x-μ)^2)) + δ*sqrt(α^2-beta^2) + β*(x-μ)
end
