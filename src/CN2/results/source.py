
import numpy as np
from scipy.signal import welch

def load_data(filepath, variable, latitude, longitude, level):
    """ Carregar dados de uma variável específica do ERA5 para uma localização e nível específico. """
    ds = xr.open_dataset(filepath)
    data = ds[variable].sel(latitude=latitude, longitude=longitude, level=level)
    return data

def calculate_psd(time_series, sampling_frequency=1):
    """ Calcular a densidade espectral de potência de uma série temporal. """
    frequencies, psd_values = welch(time_series, fs=sampling_frequency, nperseg=1024)
    return frequencies, psd_values

def identify_frequencies(psd_values, frequencies):
    """ Identificar f_L, f_l, f_t com base na análise da PSD. """
    # Este é um exemplo simplificado. Você deve adaptar este método para seu caso específico.
    f_L = frequencies[np.argmax(psd_values)]  # Exemplo: pico mais baixo de frequência
    f_l = frequencies[-1]  # Frequência mais alta disponível
    f_t = frequencies[len(frequencies) // 2]  # Frequência média
    return f_L, f_l, f_t

def main():
    # Configurar parâmetros
    filepath = 'path_to_your_era5_data.nc'  # Caminho para o arquivo NetCDF
    variable = 'u'  # Velocidade do vento na direção leste (u-component)
    latitude = 50.0  # Exemplo de latitude
    longitude = 10.0  # Exemplo de longitude
    level = 850  # Exemplo de nível de pressão em hPa

    # Carregar dados
    wind_speed = load_data(filepath, variable, latitude, longitude, level)

    # Calcular PSD
    frequencies, psd_values = calculate_psd(wind_speed.values)

    # Identificar frequências
    f_L, f_l, f_t = identify_frequencies(psd_values, frequencies)

    print(f"f_L: {f_L} Hz, f_l: {f_l} Hz, f_t: {f_t} Hz")

if __name__ == '__main__':
    main()



OTROOOOO

import numpy as np
import xarray as xr
from scipy.signal import welch

# Carregar dados ERA5 (temperatura e velocidade do vento)
ds = xr.open_dataset('path_to_ERA5_data.nc')
temperature = ds['t']  # Variável de temperatura
wind_speed = ds['u']  # Variável de velocidade do vento (u-component)

# Função para calcular a Transformada de Fourier e o espectro de potência
def calculate_power_spectral_density(data, fs):
    freqs, psd = welch(data, fs=fs, nperseg=1024)
    return freqs, psd

# Calcular o espectro de potência para a temperatura e velocidade do vento
fs = 1  # Frequência de amostragem (por exemplo, 1 medição por hora)
freqs_temp, psd_temp = calculate_power_spectral_density(temperature.values.flatten(), fs)
freqs_wind, psd_wind = calculate_power_spectral_density(wind_speed.values.flatten(), fs)

# Identificar fL, ft, fl a partir do espectro de potência
def identify_frequencies(freqs, psd):
    fL = freqs[np.argmax(psd)]  # Frequência com maior potência (baixa frequência)
    ft = freqs[np.argmin(np.abs(psd - np.mean(psd)))]  # Frequência de transição
    fl = freqs[-1]  # Maior frequência (alta frequência)
    return fL, ft, fl

fL_temp, ft_temp, fl_temp = identify_frequencies(freqs_temp, psd_temp)
fL_wind, ft_wind, fl_wind = identify_frequencies(freqs_wind, psd_wind)

# Calcular E(fL) integrando o espectro de potência na baixa frequência
def calculate_E_fL(freqs, psd, fL):
    return np.trapz(psd[freqs <= fL], freqs[freqs <= fL])

E_fL_temp = calculate_E_fL(freqs_temp, psd_temp, fL_temp)

# Mostrar os resultados
print(f"Frequências identificadas (Temperatura): fL={fL_temp}, ft={ft_temp}, fl={fl_temp}")
print(f"Frequências identificadas (Velocidade do Vento): fL={fL_wind}, ft={ft_wind}, fl={fl_wind}")
print(f"Densidade espectral de potência E(fL) (Temperatura): {E_fL_temp}")

# Implementar a fórmula completa para calcular C_n^2 usando os valores calculados
A = 80e-60  # em hPa
pressure = ds['p']  # Variável de pressão
T_mean = temperature.mean(dim='time')

C_n2 = (A * pressure / T_mean**2)**2 * E_fL_temp * np.exp(-3 * np.log(ft_temp / fL_temp) - 5/3 * np.log(fl_temp / ft_temp))**0.125 * fl_temp**(5/3)
print(C_n2)


