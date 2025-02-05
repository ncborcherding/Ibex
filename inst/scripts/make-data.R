#########################
#Defining Hyperparameters
#########################

factors <- c("OHE", "atchleyFactors", "crucianiProperties", "kideraFactors", "MSWHIM", "tScales")
hidden_dim1 <- 512
hidden_dim2 <- 256
latent_dim <- 128
batch_size <- 128
learning_rate <- 1e-6
epochs <- 128
optimizer <- "adam"
layer_act <- "relu"
epsilon.std <- 1
amino.acids <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

####################
#Training CNN Models
####################

set.seed(42)

for(i in seq_along(factors)) {
  
  es = callback_early_stopping(
    monitor = "val_loss",
    min_delta = 0,
    patience = 8,
    verbose = 1,
    mode = "min")
  
  sequence.matrix <- readRDS(paste0(data.path, factors[i], "_Heavy_CDR3.rds"))
  
  stratified.sequences <- prepare_data(sequence.matrix, 
                                       train_split = 0.75, 
                                       val_split = 0.2)
  
  # Create the training, validation, and test sets
  x_train <- stratified.sequences[[1]]
  x_val <- stratified.sequences[[2]]
  x_test <- stratified.sequences[[3]]
  rm(stratified.sequences)
  rm(sequence.matrix)
  gc()
  
  d.1 <- dim(x_train)[2]
  input_layer <- layer_input(shape = c(d.1))
  # Encoder part
  encoder <- input_layer %>%
    layer_dense(units = hidden_dim1, name = "e.1") %>%
    layer_batch_normalization(name = "bn.1") %>%
    layer_activation(activation = layer_act, name = "act.1") %>%
    layer_dense(units = hidden_dim2, name = "e.2") %>%
    layer_batch_normalization(name = "bn.2") %>%
    layer_activation(activation = layer_act, name = "act.2") %>%
    layer_dense(units = latent_dim, activation = layer_act, name = "latent_space")
  
  # Decoder part
  decoder <- encoder %>%
    layer_dense(units = hidden_dim2, name = "d.1") %>%
    layer_batch_normalization(name = "bn.3") %>%
    layer_activation(activation = layer_act, name = "act.3") %>%
    layer_dense(units = hidden_dim1, name = "d.2") %>%
    layer_batch_normalization(name = "bn.4") %>%
    layer_activation(activation = layer_act, name = "act.4") %>%
    layer_dense(units = d.1, activation = 'sigmoid', name = "output")
  
  # Complete autoencoder model
  autoencoder <- keras_model(input_layer, decoder)
  
  # Extract the latent space output
  encoder_model <- keras_model(inputs = autoencoder$input, outputs = get_layer(autoencoder, "latent_space")$output)
  
  # Create the decoder model
  latent_input <- layer_input(shape = latent_dim, name = "latent_input")
  decoder_output <- latent_input %>%
    get_layer(autoencoder, "d.1")(.) %>%
    get_layer(autoencoder, "bn.3")(.) %>%
    get_layer(autoencoder, "act.3")(.) %>%
    get_layer(autoencoder, "d.2")(.) %>%
    get_layer(autoencoder, "bn.4")(.) %>%
    get_layer(autoencoder, "act.4")(.) %>%
    get_layer(autoencoder, "output")(.)
  
  decoder_model <- keras_model(latent_input, decoder_output)
  
  autoencoder %>% compile(
    optimizer = optimizer_adam(learning_rate = learning_rate),
    loss = "mean_squared_error",
    metrics = 'mean_absolute_error')
  
  # Train the model
  history <- autoencoder %>% fit(
    x = x_train,
    y = x_train,
    validation_data = list(x_val, x_val),
    epochs = epochs,
    batch_size = batch_size,
    shuffle = TRUE,
    callbacks = es)
  
  save_model(encoder_model, paste0(data.path, "/models/Human_Heavy_CNN_", factors[i], "_encoder.keras"), overwrite = TRUE)
  save_model(decoder_model, paste0(data.path, "/models/Human_Heavy_CNN_", factors[i], "_decoder.keras"), overwrite = TRUE)
  save_model(autoencoder, paste0(data.path, "/models/Human_Heavy_CNN_", factors[i], "_autoencoder.keras"), overwrite = TRUE)
}

####################
#Training VAE Models
####################

for(i in seq_along(factors)) {
  
  es = callback_early_stopping(
    monitor = "val_loss",
    min_delta = 0,
    patience = 8,
    verbose = 1,
    mode = "min")
  
  sequence.matrix <- readRDS(paste0(data.path, factors[i], "_Heavy_CDR3.rds"))
  
  stratified.sequences <- prepare_data(sequence.matrix, 
                                       train_split = 0.75, 
                                       val_split = 0.2)
  
  # Create the training, validation, and test sets
  x_train <- stratified.sequences[[1]]
  x_val <- stratified.sequences[[2]]
  x_test <- stratified.sequences[[3]]
  rm(stratified.sequences)
  rm(sequence.matrix)
  gc()
  
  vae_loss_layer <- function(original_dim) {
    layer_lambda(
      f = function(x) {
        x_decoded_mean <- x[[1]]
        x_input <- x[[2]]
        z_mean <- x[[3]]
        z_log_var <- x[[4]]
        
        # Reconstruction loss
        xent_loss <- loss_mean_squared_error(x_input, x_decoded_mean) * original_dim
        
        # KL Divergence loss
        kl_loss <- -0.5 * tf$reduce_mean(1 + z_log_var - tf$square(z_mean) - tf$exp(z_log_var), axis = -1L)
        
        # Total loss
        tf$reduce_mean(xent_loss + kl_loss)
      },
      output_shape = list(NULL, 1)  # Explicit output shape
    )
  }
  original_dim <- ncol(x_test)
  
  
  # Encoder
  encoder_input <- layer_input(shape = original_dim)
  h <- encoder_input
  h <- layer_dense(h,
                   units = hidden_dim1, 
                   activation = layer_act, 
                   name = "e.1")
  h <- layer_dense(h,
                   units = hidden_dim2, 
                   activation = layer_act, 
                   name = "e.2")
  z_mean <- layer_dense(h, units = latent_dim, name = "z_mean")
  z_log_var <- layer_dense(h, units = latent_dim, name = "z_log_var")
  
  # Sampling Layer
  z <- layer_lambda(f = function(args) {
    z_mean <- args[[1]]
    z_log_var <- args[[2]]
    batch <- tf$shape(z_mean)[1]
    dim <- tf$shape(z_mean)[2]
    epsilon <- tf$random$normal(shape = c(batch, dim), mean = 0., stddev = epsilon.std)
    z_mean + tf$exp(z_log_var / 2) * epsilon
  }, output_shape = c(latent_dim))(list(z_mean, z_log_var))
  
  # Decoder
  decoder_input <- layer_input(shape = latent_dim)
  d <- decoder_input
  d <- layer_dense(d,
                   units = hidden_dim2, 
                   activation = layer_act, 
                   name = "d.1")
  d <- layer_dense(d,
                   units = hidden_dim1, 
                   activation = layer_act, 
                   name = "d.2")
  decoder_output <- layer_dense(d, units = original_dim, activation = "sigmoid")
  
  # Encoder and Decoder Models
  encoder <- keras_model(encoder_input, z_mean)
  decoder <- keras_model(decoder_input, decoder_output)
  
  # VAE Model
  decoder_output <- decoder(z)
  vae <- keras_model(encoder_input, decoder_output)
  
  # Add custom loss layer
  loss_layer <- vae_loss_layer(original_dim)(list(decoder_output, encoder_input, z_mean, z_log_var))
  vae_with_loss <- keras_model(encoder_input, loss_layer)
  
  # Dummy loss function
  dummy_loss <- function(y_true, y_pred) {
    tf$reduce_mean(y_pred)
  }
  
  # Compile the model
  vae_with_loss %>% compile(optimizer = optimizer_adam(learning_rate = learning_rate), 
                            loss = dummy_loss, 
                            metrics = c("mean_squared_error", "mean_absolute_error"))
  
  history <- vae_with_loss %>% fit(
    x_train, x_train, 
    shuffle = TRUE,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(x_test, x_test),
    verbose = 0,
    callbacks = es
  )
  
  save_model(encoder, paste0(data.path, "/models/Human_Heavy_VAE_", factors[i], "_encoder.keras"), overwrite = TRUE)
  save_model(decoder, paste0(data.path, "models/Human_Heavy_VAE_", factors[i], "_decoder.keras"), overwrite = TRUE)
  save_model(vae, paste0(data.path, "models/Human_Heavy_VAE_", factors[i], "_autoencoder.keras"), overwrite = TRUE)
}
