@test norm(av_a - α)/abs(α) < 0.01
@test abs(n̄ - abs(α).^2)/abs(α)^2 < 0.01
@test abs(Vn - abs(α)^2)/abs(α)^2 < 0.05
